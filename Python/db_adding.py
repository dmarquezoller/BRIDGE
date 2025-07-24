import sqlite3
import pandas as pd
import os

############# REQUIREMENTS OF YOUR DATA #############
#                                                   #
#   1. It must be a csv.                            #
#   2. You require AT LEAST 3 identifier columns:   #
#        - gene name -> renamed -> Gene_Name        #
#        - gene id -> renamed -> Gene_ID            #
#        - protein id -> renamed -> Protein ID      #
#      All of the name rules specified must be      #
#      strictly followed.                           #
#   3. All the timeline columns must have a         #
#      single integer in the end specifying         #
#      the replicate                                #
#        -i.e. X6_hpf_1                             #
#   4. There must be no NAs in any identifier       #
#      columns                                      #
#   5. The name of the table must adhere to the     #
#      following structure:                         #
#        - species_datatype_(optional info)_id       #
#      i.e. zebrafish_rnaseq_investigatorA_1        #                          
#   6. For phosphoproteomics an additional          #
#      identifier is needed: the peptide with the   #
#      mutation, and it has to be called pepG       #
#      i.e. AAAGDEAGGsSR_p1_ac0                     #                
#                                                   #
#                                                   #
#   * If any of these rules is not met for a        #
#     table submited to the database, the app       #
#     will most likely crash                        #                            
#                                                   #
#####################################################


# === Function to enable column selection ===

def prompt_column_selection(df, prompt_text):
    print(f"\n{prompt_text}")
    print("Available columns:")
    for idx, col in enumerate(df.columns):
        print(f"{idx}: {col}")
    indices = input(f"Enter {prompt_text} numbers sepparated by commas or ranges (e.g. 1,2,3,5:10): ")
    selected = []

    for part in indices.split(","):
        part = part.strip()
        if ":" in part:
            start, end = part.split(":")
            if start.strip().isdigit() and end.strip().isdigit():
                selected.extend(df.columns[int(i)] for i in range(int(start), int(end) + 1))
        elif part.isdigit():
            selected.append(df.columns[int(part)])

    return selected


def main():
    # === Step 1: Input CSV path ===
    csv_path = input("Enter path to the CSV file: ").strip()
    if not os.path.exists(csv_path):
        print("File does not exist.")
        return

    # === Step 2: Read CSV ===
    df = pd.read_csv(csv_path)
    print(f"CSV loaded with {df.shape[0]} rows and {df.shape[1]} columns.")

    # === Step 3: Ask user to select identifier and timepoint columns ===
    id_cols = prompt_column_selection(df, "identifier columns")
    tp_cols = prompt_column_selection(df, "timepoint columns")

    # === Step 4: Ask for table name ===
    table_name = input("\nEnter a name for the new table in the database: ").strip()
    if not table_name.isidentifier():
        print("Invalid table name.")
        return

    # === Step 5: Connect to SQLite ===
    db_path = input("\nEnter path to the database: ").strip() 
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # === Step 6: Upload table ===
    try:
        df.to_sql(table_name, conn, if_exists='replace', index=False)
        print(f"Table '{table_name}' uploaded to the database.")
    except Exception as e:
        print(f"Error uploading the table to the database: {e}")
        conn.close()
        return

    # === Step 7: Insert metadata ===
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS table_metadata (
            table_name TEXT PRIMARY KEY,
            identifier_columns TEXT,
            timepoint_columns TEXT
        )
    """)

    # Convert column lists to comma-separated strings
    id_str = ",".join(id_cols)
    tp_str = ",".join(tp_cols)

    cursor.execute("""
        INSERT OR REPLACE INTO table_metadata (table_name, identifier_columns, timepoint_columns)
        VALUES (?, ?, ?)
    """, (table_name, id_str, tp_str))

    conn.commit()
    conn.close()
    print("Metadata added to 'table_metadata'.")

if __name__ == "__main__":
    main()
