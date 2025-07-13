import os
import csv

def listar_csv_info(carpeta):
    """
    Recorre todos los archivos CSV en la carpeta indicada y devuelve una lista de diccionarios
    con: nombredecsv, primera línea, primeras 3 líneas como ejemplo.

    Parámetros:
    carpeta (str): Ruta de la carpeta que contiene los archivos CSV.

    Retorna:
    List[Dict[str, Any]]: Lista de diccionarios con la información solicitada.
    """
    resultados = []
    # Recorre todos los ficheros con extensión .csv
    for nombre_archivo in os.listdir(carpeta):
        if nombre_archivo.lower().endswith('.csv'):
            ruta = os.path.join(carpeta, nombre_archivo)
            try:
                with open(ruta, newline='', encoding='utf-8') as csvfile:
                    reader = csv.reader(csvfile)
                    # Leer las líneas
                    lineas = list(reader)
                    if not lineas:
                        # Archivo vacío
                        primera_linea = None
                        ejemplo_3 = []
                    else:
                        primera_linea = ','.join(lineas[0])
                        # Obtener hasta las primeras 3 líneas
                        ejemplo_3 = [','.join(fila) for fila in lineas[:3]]

                resultados.append({
                    'nombredecsv': nombre_archivo,
                    'primera_línea': primera_linea,
                    'primeras_3_líneas_como_ejemplo': ejemplo_3
                })
            except Exception as e:
                print(f"Error al procesar {ruta}: {e}")
    return resultados


def guardar_info_txt(info_csv, archivo_salida):
    """
    Guarda la información de los CSV en un archivo de texto.

    Parámetros:
    info_csv (list): Lista de diccionarios con la información de los CSV.
    archivo_salida (str): Ruta del archivo de texto de salida.
    """
    with open(archivo_salida, 'w', encoding='utf-8') as f:
        for info in info_csv:
            f.write(f"Nombre del CSV: {info['nombredecsv']}\n")
            f.write(f"Primera línea: {info['primera_línea']}\n")
            f.write("Primeras 3 líneas como ejemplo:\n")
            for linea in info['primeras_3_líneas_como_ejemplo']:
                f.write(f"  - {linea}\n")
            f.write("\n")

if __name__ == '__main__':
    carpeta_datos = 'Datos'
    archivo_salida = 'resumen_csv.txt'

    info_csv = listar_csv_info(carpeta_datos)
    guardar_info_txt(info_csv, archivo_salida)
    print(f"Información guardada en {archivo_salida}")
