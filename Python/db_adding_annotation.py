import sqlite3
import pandas as pd
import os

############# REQUIREMENTS OF YOUR DATA #############
#                                                   #
#   1. It must be a csv.                            #
#   2. You require the following columns:           #
#          - Gene_ID                                #
#          - Gene_Name                              #
#          - Chromosome                             #
#          - Gene_Start                             #
#          - Gene_End                               #
#          - Gene_Type                              #
#          - Strand                                 #
#   3. The names of the columns have to be          #
#      the precise ones above.                      #
#   4. You can download this annotation table       #
#      with these columns in BioMart (ENSEMBL).     #          #
#   5. The name of the annotation table uploaded    #
#      has to follow the strict rule:               #
#      <specie>_annotation                          #
#      i.e. zebrafish_annotation                    #
#   6. The name of the specie has to be exactly     #
#      the same as the one used for the data        #
#      table.                                       #
#                                                   #
#                                                   #
#   * If any of these rules is not met for a        #
#     table submited to the database, the app       #
#     will most likely crash                        #                            
#                                                   #
#####################################################


def load_csv_with_validation():
    file_path = input("Enter the path to the annotation CSV file: ").strip()
    if not os.path.exists(file_path):
        print("File does not exist.")
        return None, None

    try:
        df = pd.read_csv(file_path)
        print(f"Annotation CSV loaded with {df.shape[0]} rows and {df.shape[1]} columns.")
        return df, file_path
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return None, None

def upload_dataframe_to_sql(df, conn, table_name):
    try:
        df.to_sql(table_name, conn, if_exists='replace', index=False)
        print(f"Annotation table '{table_name}' uploaded to the database.")
        return True
    except Exception as e:
        print(f"Error uploading the table to the database: {e}")
        return False

def main():
    print("=== Annotation Upload Script ===")

    # === Load annotation CSV ===
    df, file_path = load_csv_with_validation()
    if df is None:
        return

    # === Ask for table name ===
    table_name = input("Enter a name for the annotation table in the database: ").strip()
    if not table_name.isidentifier():
        print("Invalid table name (must be alphanumeric with no spaces or symbols except underscores).")
        return

    # === Connect to database ===
    db_path = input("\nEnter path to the database: ").strip()  
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # === Upload annotation table ===
    if not upload_dataframe_to_sql(df, conn, table_name):
        conn.close()
        return

    # === Store annotation metadata ===
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS annotation_metadata (
            table_name TEXT PRIMARY KEY,
            source_file TEXT,
            row_count INTEGER,
            column_count INTEGER
        )
    """)
    cursor.execute("""
        INSERT OR REPLACE INTO annotation_metadata (table_name, source_file, row_count, column_count)
        VALUES (?, ?, ?, ?)
    """, (table_name, os.path.basename(file_path), df.shape[0], df.shape[1]))
    conn.commit()

    conn.close()
    print("Annotation metadata saved to 'annotation_metadata'. Done.")

if __name__ == "__main__":
    main()
