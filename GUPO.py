from tkinter import *
from tkinter.filedialog import askopenfilename,askdirectory
import os
import shutil
import time
import csv
import wget
from operator import itemgetter
import requests
import json
import magic
import zipfile

errores=[] # Para almacenar los errores que se vayan produciendo
experimentos=[] # Para los experimentos seleccionados realmente. Tuplas (exp acc, comparison)
archivos=[] # Para los archivos descargados realmente (nombre - path). El nombre coincide con el exp acc
gen_principal="" # El gen objeto de estudio
DESCONOCIDO="Desconocido" # "Constante" para cuando no podemos saber qué gen principal estamos tratando
resultado={} # Diccionario en el que la key es un cuarteto ((id,gen),(exp,comp)) y el valor es una tupla ("valor de expresión",pvalue)
resultado_gen_ppal={} # Lo mismo, pero sólo con los valores del gen principal y clave un par (exp,comp)
genes_secundarios=[] # Para los genes resultantes. El primer ídice de "resultado"
parejas=[] # Para los pares exp-comparison. Segundo índice de "resultado"
dir_trabajo=""
dir_ejec = ""
dir_datos = ""
dir_exp = ""


########################################################################################################################
#                                                                                                                      #
#                                               FUNCIONES                                                              #
#                                                                                                                      #
########################################################################################################################

################################# Muestra diálogo de selección de archivo ##############################################
def busca_archivo():
    archivo = askopenfilename()
    if archivo!="":
        fuente_texto.delete(0, END)
        fuente_texto.insert(0, archivo)

################################# Muestra diálogo de selección de directorio############################################
def busca_directorio():
    directorio=askdirectory()
    if directorio!="":
        directorio_texto.delete(0,END)
        directorio_texto.insert(0,directorio)

################################# Envía un texto al log de ejecución ###################################################
def muestra_log(texto,seguido=FALSE):
    log.config(state=NORMAL)
    if seguido:
        log.insert(END, texto + " ")
    else:
        log.insert(END, texto+"\n")
    log.config(state=DISABLED)

    log.see(END) # Para ver siempre lo último mostrado

    root.update_idletasks() # FUNDAMENTAL !!!!!!!!!  para que refresque el GUI sin esperar a que se haga todo

############################# Guarda un error en la lista y lo envía al log ############################################
def registra_error(texto):
    errores.append(texto)
    muestra_log(texto)

######################### Construye un texto precedido de la hora (con línea en blanco) para el log ####################
def ahora(texto):
    return ("\n"+time.strftime("%H:%M:%S")+" - "+texto)

################## Muestra el mensaje de finalización de ejecución y habilita el botón de ejecutar #####################
def finaliza():
    muestra_log(ahora("Ejecución Finalizada"))
    ejecutar_boton.config(state=NORMAL)
    fuente_boton.config(state=NORMAL)
    directorio_boton.config(state=NORMAL)

############################# Pasos previos y posteriores al proceso de los datos ######################################
def ejecutar():
    # Vaciamos las listas que vamos a usar
    errores.clear()
    experimentos.clear()
    archivos.clear()
    resultado.clear()
    resultado_gen_ppal.clear()
    genes_secundarios.clear()
    parejas.clear()
    # Deshabilita el botón ejecutar para que no se pueda pulsar durante la ejecución
    ejecutar_boton.config(state=DISABLED)
    fuente_boton.config(state=DISABLED)
    directorio_boton.config(state=DISABLED)
    # limpia el log de ejecución
    log.config(state=NORMAL)
    log.delete(1.0,END)
    log.config(state=DISABLED)
    # Si hay inconsistencias en los datos de la pantalla de selección, cancelo ejecución
    if errores_previos():
        ejecutar_boton.config(state=NORMAL)
        fuente_boton.config(state=NORMAL)
        directorio_boton.config(state=NORMAL)
        return
    # Si todo ha ido bien, proceso los datos ...
    procesar()
    # Genero los archivos de Resultado y Errores
    genera_archivos()
    # Y se acabó
    finaliza()

######################### Chequea que los datos de la pantalla de selección sean coherentes ############################
def errores_previos():
    global gen_principal
    url_1="http://www.ebi.ac.uk/ebisearch/ws/rest/atlas-genes?query=" #Para la llamad al API
    url_2="&fields=name&format=json"
    ids = []
    genes = []
    # El campo de archivo está informado?
    if fuente_texto.get()=="":
        muestra_log("  Especifica un archivo")
        return TRUE
    # El archivo existe?
    elif os.path.isfile(fuente_texto.get())==FALSE:
        muestra_log("  El archivo indicado no existe")
        return TRUE
    # Chequeos sobre la estructura y contenido del archivo
    else:
        with open(fuente_texto.get()) as fuente:
            lector = csv.reader(fuente, delimiter="\t")
            num_linea=0
            for linea in lector:
                num_linea+=1
                if num_linea==1: #La primera línea del archivo tiene 7 campos separados por tabuladores?
                    if len(linea)!=7:
                        muestra_log("  El archivo no tiene la estructura adecuada")
                        return TRUE
                else: # De cuántos genes es el archivo de búsqueda? (recojo los id para obtener los genes)
                    if linea[0] not in ids:
                        ids.append(linea[0])
        try: # Si fallo en la recuperación del gen principal lo marco mo desconocido
            for id in ids: #Para cada id recupero su gen con una llamada al API de EBI (https://www.ebi.ac.uk/ebisearch/swagger.ebi)
                url_completa=url_1+id+url_2
                respuesta=requests.get(url_completa)
                mi_json=json.loads(respuesta.text)
                mi_lista= mi_json["entries"]
                for l in mi_lista: # De las diferentes entradas me quedo la que coincide el id y cojo su gen (en mays)
                    if l["id"]==id:
                        gen=str(l["fields"]["name"][0]).upper()
                        if gen not in genes:
                            genes.append(gen)
            if len(genes)>1: # Sólo se aceptan archivos de 1 gen
                muestra_log("  Error en archivo de búsqueda de Expression Atlas.\n"
                            "  Debería corresponder a un único gen. Se han encontrado "+str(len(genes))+":")
                for g in genes:
                    muestra_log("  "+g)
                return TRUE
            else:
                gen_principal=genes[0]
        except:
            gen_principal=DESCONOCIDO
    # El campo de directorio está informado?
    if directorio_texto.get()=="":
        muestra_log("  Especifica un directorio")
        return TRUE
    elif os.path.isdir(directorio_texto.get())==FALSE: # El directorio existe?
        muestra_log("  El directorio indicado no existe")
        return TRUE
    # El valor de expresión introducido es numérico?
    try:
        float(exp_valor.get())
    except:
        muestra_log ("  El umbral de expresión debe ser numérico")
        return TRUE
    if float(exp_valor.get())<1:
        muestra_log("  El umbral de expresión debe ser >=1")
        return TRUE
    # El valor de pValue máximo está en [0,1]?
    try:
        pm=float(pval_valor.get())
    except:
        muestra_log ("  El pValue máximo debe ser numérico")
        return TRUE
    if pm>1 or pm<0:
        muestra_log("  El pValue máximo debe estar entre 0 y 1")
        return TRUE
    elif pm!=0 and pm<1e-307:
        muestra_log("  El mínimo valor de pValue, distinto de cero, permitido es 1e-307")
        return TRUE

    return FALSE

############################# Ejecución secuencial de los diferentes pasos del proceso #################################
def procesar():
    muestra_log(ahora("Inicio de ejecución ..."))
    if gen_principal==DESCONOCIDO:
        muestra_log(" Gen objeto de estudio: "+gen_principal+" Error de enlace con API de EBI")
    else:
        muestra_log("  Gen objeto de estudio: " + gen_principal)
    muestra_log("  Umbral de expresión: "+str(abs(float(exp_valor.get()))))
    muestra_log("  pVal máximo: "+str(float(pval_valor.get())))
    muestra_log("  Especie: " + vpd.get())
    muestra_log("  Generación de archivos adicionales: "+vfaa.get())
    # Primero, crear los directorios y copiar a ellos el archivo fuente
    if crea_directorios() == FALSE:
        return
    # Segundo, filtrar los experimentos que hay en el archivo fuente
    filtra_fuente()
    if (len(experimentos) == 0):
        muestra_log("  No hay experimentos que cumplan las condiciones expecificadas")
        return
    # Tercero, descargar los experimentos que han pasado el filtro
    descarga_experimentos()
    if (len(archivos)==0):
        muestra_log("  No se podido descargar ningún archivo")
        return
    # Finalmente, revisa los experimentos para obtener la lista de genes sobre/sub expresados
    revisa_experimentos()
    return

########################## Crea el directorio de ejecución dentro del de trabajo #######################################
def crea_directorios():
    global dir_trabajo,dir_ejec,dir_datos,dir_exp
    # El nombre será: fechahora formato yankee con el nombre del gen, el umbral, el pvalue, y el flag de archivo adicional
    umbral=str(abs(float(exp_valor.get())))
    pvalue=str(float(pval_valor.get()))
    dir_trabajo=directorio_texto.get()
    # dir_ejec = dir_trabajo+"/"+time.strftime("%Y%m%d") + time.strftime("%H%M%S")+"_"+gen_principal
    dir_ejec=dir_trabajo+"/"+gen_principal+"("+umbral+"_"+pvalue+"_"+vpd.get()+"_"+vfaa.get()+")"+"_"+time.strftime("%Y%m%d") + time.strftime("%H%M%S")
    dir_datos = dir_ejec + "/ORIGINAL"
    dir_exp = dir_ejec + "/EXPERIMENTOS"
    try:
        os.makedirs(dir_ejec) # Directorio de ejecución y dos subdirectorios
        os.makedirs(dir_datos)
        os.makedirs(dir_exp)
        shutil.copyfile(fuente_texto.get(),dir_datos+"/"+gen_principal+".tsv") # Copia el archivo original de Expression Atlas
        muestra_log("  Creado directorio de ejecución: "+dir_ejec)
        return TRUE
    except:
        muestra_log("  Error al crear el directorio de ejecución")
        return FALSE

############################# Filtra los experimentos que hay en el archivo fuente #####################################
def filtra_fuente():
    muestra_log(ahora("Filtrando experimentos de búsqueda ..."))
    root.update_idletasks()
    num_linea=0
    p_valor=abs(float(exp_valor.get())) #valor de expresión introducido en pantalla
    normalizado=float(pval_valor.get())
    p_especie=especie_lista.get() #especie introducida por pantalla
    with open(fuente_texto.get()) as fuente:
        lector = csv.reader(fuente, delimiter="\t") # Los datos de la búsqueda de AtlasExpress vienen separados por comas
        for linea in lector:
            num_linea+=1
            if num_linea== 1: # Me salto la primera línea porque contiene las cabeceras
                continue
            try:
                especie=linea[1] # la especie que figura en el archivo
                foldchange=float(linea[4])
                pvalue=float(linea[5])
                if (foldchange>=p_valor or foldchange<=-p_valor) and (p_especie=="TODAS" or p_especie==especie) and (pvalue<normalizado):
                    tupla=(linea[2],linea[3])
                    if (tupla not in experimentos):  # para evitar entradas repetidas
                        experimentos.append(tupla)
                        muestra_log("  Experimento seleccionado: "+tupla[0]+" - "+tupla[1])
            except:
                registra_error("  Formato erróneo en línea "+str(num_linea)+" del archivo de Expression Atlas")

########################## Descarga los experimentos que han pasado el filtro ##########################################
def descarga_experimentos():
    muestra_log(ahora("Descarga de archivos de experimentos seleccionados:"))
    patron=r'zip' # regular expression para comprobar el tipo de archivo descargado
    descargados=[] # para no repetir desacarga de experimentos (con mismo exp acc y diferente comparison)
    viejos=[] # Experimentos asociados a archivos zip
    nuevos=[] # experimentos asociados a componentes de los zip
    for exp in experimentos:
        if exp[0] in descargados:
            continue
        muestra_log("  Descargando archivo " + exp[0] + " ... ", TRUE)
        root.update_idletasks()
        destino=dir_exp+"/"+exp[0]
        resultado_descarga = descarga_individual(exp[0], destino)
        if resultado_descarga == "ok":
            descargados.append(exp[0])
            tipo_mime=magic.from_file(destino,mime=True)
            if re.search(patron,tipo_mime): # Si el archivo descargado es zip -> debo procesar su contenido
                componentes=procesa_zip(dir_exp,destino)
                if len(componentes)==0: # Si ha fallado al descomprimir, paso de él
                    muestra_log("ERROR al descomprimir!")
                    continue
                comparison = exp[1] # Si ha ido bien ...
                for c in componentes:
                    nombre_completo=dir_exp+"/"+c
                    tupla = (c, nombre_completo)
                    archivos.append(tupla) # ... considero cada componente como un archivo ...
                    nuevos.append((c,comparison)) # ... guardo los experimentos asociados ...
                if exp not in viejos:
                    viejos.append(exp) # ... y marco para borrar el del propio zip
            else:
                tupla = (exp[0], destino)
                archivos.append(tupla)
        muestra_log(resultado_descarga)
    # Una vez que termino de procesar los experimentos, sustituir los de los zip por los de sus componentes
    for exp in viejos:
        experimentos.remove(exp)
    for exp in nuevos:
        experimentos.append(exp)

###################### Pruebo a descargar una archivo de dos BBDD diferentes conocidas #################################
def descarga_individual(nombre,path_completo):
    inicio_url = "https://www.ebi.ac.uk/gxa/experiments-content/"
    fin_url_Microarray = "/resources/ExperimentDownloadSupplier.Microarray/query-results"
    fin_url_RnaSeq = "/resources/ExperimentDownloadSupplier.RnaSeqDifferential/tsv"
    try:  # Primero busco el archivo entre los experimentos de Microarray
        nombre_url = inicio_url + nombre + fin_url_Microarray
        wget.download(nombre_url, path_completo)
        return ("ok")
    except:
        try:  # Si no está,lo busco entre los de RnaSeq
            nombre_url = inicio_url + nombre + fin_url_RnaSeq
            wget.download(nombre_url, path_completo)
            return ("ok")
        except Exception as e:
            return("ERROR! -> " + str(e))

############################ Sustituye un archivo comprimido por su contenido ##########################################
def procesa_zip(directorio,path_completo):
    contenido=[]
    nombres=[]
    try:
        mi_zip=zipfile.ZipFile(path_completo,mode='r')
        contenido=mi_zip.namelist() # Averiguo los nombres de sus componentes
        mi_zip.extractall(path=directorio) # Extraigo el contenido en el mismo directorio
        for c in contenido: # Elimino la extensión del nombre de cada archivo extraido para que estén como los demás
            n=os.path.splitext(c)[0]
            original=directorio+"/"+c
            nuevo=directorio+"/"+n
            os.rename(original,nuevo)
            if n not in nombres:
                nombres.append(n)
            break # ME HAN PEDIDO QUE DE LOS ZIP, SÓLO COJA UNO DE LOS COMPONENTES !!!!!!!!!!!!!
        mi_zip.close()
        os.remove(path_completo)  # Y borro el zip
    except: # Para poder detectar los errores a la vuelta, si algo ha fallado, devuelvo una lista de archivos vacía
        nombres.clear()
    return nombres

################ Revisa los datos de los experimentos para encontrar genes sobre/sub expresados ########################
def revisa_experimentos():
    muestra_log(ahora("Tratamiento de los archivos descargados:"))
    # El objetivo es construir un array bidimensional asociativo en ambos índices:
    # Una tabla que tenga como índice de fila "el nombre del gen",
    # como índice de columna el nombre del experimento, o sea, "accesion + comparison",
    # y en su intersección una tupla: el foldchange del gen en dicho experimento (si es que aparece sobre/sub) expresado y el pvalue
    # Luego se volcará este array en un csv para abrir con hoja de cálculo, que es el objeto final del programa
    # Ese array bidim se hará con un diccionario que tiene como clave una tupla (gen,(exp acc,comp))

    for tupla in archivos:
        exp=tupla[0]
        muestra_log("  Procesando archivo " + exp + " ... ", TRUE)
        columnas=localiza_columnas(exp) # Nombres de las columnas que debo revisar (versión reducida de "experimentos",i.e., tuplas(exp,comp))
        indices={} # Posiciones de las columnas en el archivo del exper. Diccionario: clave=comp, valor=(pos foldchange, pos pvalue)
        try:
            archivo=open(tupla[1],"r")
            lector=csv.reader(archivo, delimiter="\t") # Abro el archivo para procesarlo como csv
            contador_linea = 0
            for linea in lector:
                contador_linea += 1
                if contador_linea>4: # Estas son las líneas que continen los datos del experimento para cada gen
                    gen=(linea[0],linea[1])
                    if gen not in genes_secundarios:
                        genes_secundarios.append(gen)
                    for k,v in indices.items():
                        par=(exp,k) # accesion,comparison
                        if par not in parejas:
                            parejas.append(par)
                        cuarteto=(gen,par)
                        try:
                            foldchange=float(linea[v[0]])
                        except:
                            foldchange=0.99
                        try:
                            pvalue=float(linea[v[1]])
                        except:
                            pvalue=1.01
                        # # Actualizamos la megamatriz de resultados
                        # if cuarteto in resultado:
                        #     if abs(resultado[cuarteto][0])<abs(foldchange): # si tiene varios foldchange me quedo con el de mayor abs
                        #         resultado[cuarteto]=(foldchange,pvalue)
                        # else:
                        #     resultado[cuarteto]=(foldchange,pvalue)
                        # # Actualizamos la minimatriz con los resultados específicos del gen principal (una sola línea)
                        # if cuarteto[0][1]==gen_principal:
                        #     trio=(gen_principal,cuarteto[1])
                        actualiza_resultados(cuarteto,foldchange,pvalue)

                elif contador_linea<4: # Las tres primeras líneas son de comentarios
                    continue
                else: # La cuarta línea contiene los nombres de las columnas (entre ellos, los comparison)
                    for col in columnas:
                        comparison=col[1]
                        i=-1 # posición del foldchange de ese comparison
                        j=-1 # posición del pvalue de ese comparison
                        for n in range(len(linea)):

                            lineasinespacios=linea[n].replace(" ", "")
                            comparisonmasfoldchange=comparison+".foldChange"
                            comparisonmaspvalue=comparison+".pValue"
                            if lineasinespacios==comparisonmasfoldchange.replace(" ", ""):
                                i=n
                            elif lineasinespacios==comparisonmaspvalue.replace(" ", ""):
                                j=n
                            # GAG: he tenido que corregirlo porque los títulos de las columnas, a veces, incluyen " " alrededor del "."
                            # if linea[n]==comparison+".foldChange":
                            #     i=n
                            # elif linea[n]==comparison+".pValue":
                            #     j=n
                        indices[comparison]=(i,j)
            muestra_log("ok")
        except Exception as e:
            muestra_log("ERROR! -> "+str(e))
        finally:
            archivo.close()
    genes_secundarios.sort()
    parejas.sort()  # ordenación por defecto (primero por el primer componente de la tupla y luego por el segundo)

################### Para un archivo de experimento, localiza las columnas que hay que revisar ##########################
def localiza_columnas(nombre):
    # Devuelve una lista reducida con los experimentos de ESE nombre
    columnas=[]
    for e in experimentos:
        if e[0]==nombre:
            columnas.append(e)
    return columnas

################ Actualiza las matrices de resultados: la general y la específica del gen principal ####################
def actualiza_resultados(cuarteto,foldchange,pvalue):
    # Actualizamos la megamatriz de resultados
    #muestra_log(cuarteto[0][0] + " - " + cuarteto[0][1] + " - " + cuarteto[1][0] + " - " + cuarteto[1][1] + " - "+str(foldchange)+" - "+str(pvalue), FALSE)
    if cuarteto in resultado:
        if abs(resultado[cuarteto][0])<abs(foldchange): # si tiene varios foldchange me quedo con el de mayor abs
            resultado[cuarteto]=(foldchange,pvalue)
    else:
        resultado[cuarteto]=(foldchange,pvalue)
    # Actualizamos la minimatriz con los resultados específicos del gen principal (una sola línea)
    if cuarteto[0][1]==gen_principal:
        par=cuarteto[1]
        if par in resultado_gen_ppal:
            if abs(resultado_gen_ppal[par][0])<abs(foldchange):
                resultado_gen_ppal[par]=(foldchange,pvalue)
                #muestra_log(cuarteto[0][0]+" - "+cuarteto[0][1]+" - "+cuarteto[1][0]+" - "+cuarteto[1][1],FALSE)
        else:
            resultado_gen_ppal[par]=(foldchange,pvalue)
            #muestra_log(cuarteto[0][0] + " - " + cuarteto[0][1] + " - " + cuarteto[1][0] + " - " + cuarteto[1][1],FALSE)

######################### Genera el archivo con los resultados y, si procede, otro con los errores #####################
def genera_archivos():
    muestra_log(ahora("Generación de archivos:"))
    umbral=abs(float(exp_valor.get())) # valor de expresión introducido en pantalla
    maxval=abs(float(pval_valor.get())) # valor máximo de pValue introducido en pantalla
    # Primero el archivo con el resultado
    genera_resultado_principal(umbral,maxval)
    # Luego, el archivo secundario, si así lo decide el usuario
    if vfaa.get() == "Sí":
        genera_resultado_signo(umbral,maxval)
    # Otro con el log que se ha visto por pantalla
    try:
        archivo = open(dir_ejec + "/log.txt", "w")
        archivo.write(log.get(1.0, END))
        muestra_log("  Log de actividad ... ok")
    except Exception as e:
        muestra_log("  Log de actividad ... ERROR! -> " + str(e))
    finally:
        archivo.close()


########################### Genero archivo con la matriz de resultados completa ########################################
def genera_resultado_principal(umbral,maxval):
    muestra_log("  Resultados ... ", TRUE)
    try:
        archivo= open(dir_ejec+"/resultados.tsv","w")
        wr=csv.writer(archivo,delimiter="\t")
        cabecera1=["",""] # Primera línea de cabecera
        for par in parejas:
            cabecera1.append(par[0])
        wr.writerow(cabecera1)
        cabecera2 = ["",""] # Segunda línea de cabecera
        for par in parejas:
            cabecera2.append(par[1])
        wr.writerow(cabecera2)
        for gen in genes_secundarios: # Linea de valores de expresión para cada gen
            linea=[]
            linea_valida=FALSE
            linea.append(gen[0])
            linea.append(gen[1])
            for par in parejas:
                cuarteto=(gen,par)
                if cuarteto in resultado:
                    foldchange=resultado[cuarteto][0] # resultado[cuarteto] es un par de valores(foldchange,pvalue)
                    pvalue=resultado[cuarteto][1]
                    linea.append(foldchange)
                    if abs(foldchange)>=umbral and pvalue<=maxval:
                        linea_valida=TRUE
                else:
                    linea.append("")
            if linea_valida: # La línea de ese gen sólo se mostrará si en alguno de los experimentos cumple las condiciones de foldchange y pvalue
                wr.writerow(linea)
        muestra_log("ok")
    except Exception as e:
        muestra_log("ERROR! -> " + str(e))
    finally:
        archivo.close()

########## Genera archivo con resultados que coinciden en signo (almenos en una columna) con el gen principal ##########
def genera_resultado_signo(umbral,maxval):
    muestra_log("  Resultados 2 ... ", TRUE)
    try:
        archivo = open(dir_ejec + "/resultados2.tsv", "w")
        wr = csv.writer(archivo, delimiter="\t")
        cabecera1 = ["", ""]  # Primera línea de cabecera
        for par in parejas:
            cabecera1.append(par[0])
        wr.writerow(cabecera1)
        cabecera2 = ["", ""]  # Segunda línea de cabecera
        for par in parejas:
            cabecera2.append(par[1])
        wr.writerow(cabecera2)
        for gen in genes_secundarios:  # Linea de valores de expresión para cada gen
            linea = []
            linea_valida = FALSE
            linea.append(gen[0])
            linea.append(gen[1])
            for par in parejas:
                cuarteto = (gen, par)
                if cuarteto in resultado:
                    foldchange = resultado[cuarteto][0]  # resultado[cuarteto] es un par de valores(foldchange,pvalue)
                    pvalue = resultado[cuarteto][1]
                    linea.append(foldchange)
                    if abs(foldchange) >= umbral and pvalue <= maxval: # Si supero el umbral ...
                        if foldchange*resultado_gen_ppal[par][0]>=0: # ... coincidiendo en signo con el gen principal
                            linea_valida = TRUE
                else:
                    linea.append("")
            if linea_valida:  # La línea de ese gen sólo se mostrará si en alguno de los experimentos cumple las condiciones de foldchange y pvalue
                wr.writerow(linea)
        muestra_log("ok")
    except Exception as e:
        muestra_log("ERROR! -> " + str(e))
    finally:
        archivo.close()

########################################################################################################################
#                                                                                                                      #
#                                               LAYOUT                                                                 #
#                                                                                                                      #
########################################################################################################################

root = Tk()
root.geometry("1000x600")
root.resizable(width=FALSE,height=FALSE)
root.title("Búsqueda de Genes (UPO)")

#Antes tenía tb un valor "TODAS" al principio. Lo he quitado de aquí, aunque no de la logica (por si decido volver a incluirlo)
listaespecies=["aegilops tauschii","amborella trichopoda","anas platyrhynchos","anolis carolinensis",
               "anopheles gambiae","arabidopsis lyrata","arabidopsis thaliana","aspergillus fumigatus","beta vulgaris",
               "bos taurus","brachypodium distachyon","brassica napus","brassica oleracea","brassica rapa",
               "caenorhabditis elegans","canis familiaris","chlamydomonas reinhardtii","chlorocebus sabaeus",
               "chondrus crispus","ciona intestinalis","ciona savignyi","corchorus capsularis","cyanidioschyzon merolae",
               "danio rerio","dasypus novemcinctus","drosophila melanogaster","equus caballus","galdieria sulphuraria",
               "gallus gallus","glycine max","gorilla gorilla","homo sapiens","hordeum vulgare","leersia perrieri",
               "macaca mulatta","medicago truncatula","monodelphis domestica","mus musculus","musa acuminata",
               "ornithorhynchus anatinus","oryctolagus cuniculus","oryza barthii","oryza brachyantha","oryza glaberrima",
               "oryza glumaepatula","oryza longistaminata","oryza meridionalis","oryza nivara","oryza punctata",
               "oryza rufipogon","oryza sativa","ostreococcus lucimarinus","ovis aries","pan troglodytes",
               "papio anubis","physcomitrella patens","pongo abelii","populus trichocarpa","prunus persica",
               "rattus norvegicus","saccharomyces cerevisiae","schistosoma mansoni","schizosaccharomyces pombe",
               "selaginella moellendorffii","setaria italica","solanum lycopersicum","solanum tuberosum",
               "sorghum bicolor","sus scrofa","tetraodon nigroviridis","theobroma cacao","trifolium pratense",
               "triticum aestivum","triticum urartu","vitis vinifera","xenopus tropicalis","yarrowia lipolytica",
               "zea mays"]

###################################### Sección de Parámetros ###########################################################
# Marco
parametros=LabelFrame(root,text="Parámetros de ejecución",labelanchor=NW)
parametros.place(relx=0.5,rely=0.02,anchor=N,relwidth=.99,relheight=.30)
# Selección del archivo de búsqueda descargado de AtlasExpress
fuente_etiqueta = Label(parametros, text="Archivo de búsqueda de Expression Atlas")
fuente_etiqueta.place(relx=0,rely=0.1,anchor=W)
fuente_texto=Entry(parametros,width=70)
fuente_texto.place(relx=0.3,rely=0.1,anchor=W)
fuente_boton=Button(parametros,text="Buscar",command=busca_archivo)
fuente_boton.place(relx=0.9,rely=0.1,anchor=W)
# Selección del directorio de trabajo con valor por defecto
directorio_etiqueta=Label(parametros, text="Directorio para almacenar los resultados")
directorio_etiqueta.place(relx=0,rely=.3,anchor=W)
directorio_texto=Entry(parametros,width=70)
directorio_texto.place(relx=0.3,rely=0.3,anchor=W)
directorio_texto.delete(0,END)
directorio_texto.insert(0,os.path.abspath(os.path.curdir))
directorio_boton=Button(parametros,text="Buscar",command=busca_directorio)
directorio_boton.place(relx=.9,rely=.3,anchor=W)
# Valor de sobre/sub expresión
exp_etiqueta = Label(parametros, text="Umbral de (sobre/sub) expresión deseado")
exp_etiqueta.place(relx=0,rely=0.5,anchor=W)
exp_valor=Entry(parametros,width=5)
exp_valor.place(relx=.3,rely=.5,anchor=W)
# Especie
vpd=StringVar() # valor por defecto para la lista de especies
especie_etiqueta = Label(parametros, text="Especie objeto de experimento")
especie_etiqueta.place(relx=.45, rely=0.5, anchor=W)
especie_lista = Spinbox(parametros, values=listaespecies[::-1], textvariable=vpd, wrap=TRUE)
vpd.set("homo sapiens")
especie_lista.place(relx=.65, rely=.5, anchor=W)
especie_lista.config(width=25, state="readonly")
# Valor de pValue máximo
pval_etiqueta = Label(parametros, text="pValue máximo")
pval_etiqueta.place(relx=0,rely=0.7,anchor=W)
pval_valor=Entry(parametros,width=8)
pval_valor.place(relx=.3,rely=.7,anchor=W)
pval_valor.delete(0,END)
pval_valor.insert(0,0.00001)
# Checkbox para la generación de archivos adicionales
vfaa=StringVar() # valor flag archivos adicionales
cb=Checkbutton(parametros,text=" Generar archivos adicionales",variable=vfaa,onvalue="Sí",offvalue="No")
vfaa.set("Sí")
cb.place(relx=.445,rely=.7,anchor=W)
# Botón de ejecución
ejecutar_boton=Button(parametros,text="Ejecutar",command=ejecutar) #lanza la función "ejecutar"
ejecutar_boton.place(relx=.5,rely=.9,anchor=CENTER,width=900)

###################################### Sección de Salida ###############################################################
# Marco
salida=LabelFrame(root,text="Resultado de la ejecución",labelanchor=NW)
salida.place(relx=0.5,rely=0.35,anchor=N,relwidth=.99,relheight=.63)
# Barra de scroll
deslizable = Scrollbar(salida)
deslizable.place(relx=.98, rely=.5, relwidth=.03, relheight=.95, anchor=CENTER)
# Texto
log=Text(salida,bg="#000000",fg="#00FF00",font=("Courier",9))
log.place(relx=0.48,rely=0.5,relwidth=0.95, relheight=.95,anchor=CENTER)
log.config(wrap=WORD)
# Conexión texto-barra
log.config(yscrollcommand=deslizable.set)
deslizable.config(command=log.yview)
log.config(state=DISABLED)

###################################### NO OLVIDAR! es fundamental ######################################################
root.mainloop()

