import argparse
import gzip
import os
import sqlite3

import pandas as pd

############# REQUIREMENTS OF YOUR DATA #############
#                                                   #
#   1. It must be a csv.                            #
#   2. You require AT LEAST 3 identifier columns:   #
#        - gene name -> renamed -> Gene_Name        #
#        - gene id -> renamed -> Gene_ID            #
#        - protein id -> renamed -> Protein ID      #
#      All of the name rules specified must be      #
#      strictly followed.                           #
#   3. All the value columns must have an           #
#      integer in the end specifying                #
#      the replicate, preceded by a "_".            #
#      No more underscores must be used in the name.#
#      For more separations, use other symbols.     #
#        -i.e. X6.hpf_1                             #
#   4. There can not be NAs in any identifier       #
#      columns                                      #
#   5. The name of the table must adhere to the     #
#      following structure:                         #
#        - species_datatype_(optional info)_id      #
#      i.e. zebrafish_rnaseq_investigatorA_1        #
#   6. For phosphoproteomics an additional          #
#      identifier is needed: the peptide with the   #
#      mutation, and it has to be called pepG       #
#      i.e. AAAGDEAGGsSR_p1_ac0                     #
#   7. If adding your own processed data ensure     #
#      that its an object of class                  #
#      SummarizedExperiment and that it has         #
#      been done with the same columns that the     #
#      raw table has (for matching reasons          #
#      between cache key and tables).               #
#                                                   #
#                                                   #
#   * If any of these rules is not met for a        #
#     table submited to the database, the app       #
#     will most likely crash. Read them carefuly    #
#     and ensure that your database meets all the   #
#     requirements                                  #
#                                                   #
#   IMPORTANT                                       #
#      R packages required:                         #
#       - storr, DBI, RSQLite                       #
#                                                   #
#####################################################


# === Function to enable column selection ===


def prompt_column_selection(df, prompt_text):
    print(f"\n{prompt_text}")
    print("Available columns:")
    for idx, col in enumerate(df.columns):
        print(f"{idx+1}: {col}")
    indices = input(
        f"Enter {prompt_text} numbers separated by commas or ranges using ':'. (e.g. 1,2,3,5:10): "
    )
    selected = []

    for part in indices.split(","):
        part = part.strip()
        if ":" in part:
            start, end = part.split(":")
            if start.strip().isdigit() and end.strip().isdigit():
                selected.extend(
                    df.columns[int(i) - 1] for i in range(int(start), int(end) + 1)
                )
        elif part.isdigit():
            selected.append(df.columns[int(part) - 1])

    return selected


def main():
    parser = argparse.ArgumentParser(description="Add a table to the BRIDGE database.")
    parser.add_argument("--csv", help="Path to the CSV/TSV file (can be .gz)")
    parser.add_argument("--db", help="Path to the SQLite database")
    parser.add_argument("--table", help="Name for the new table in the database")
    parser.add_argument(
        "--id-cols", help="Comma-separated indices for identifier columns (e.g. 0,1,2)"
    )
    parser.add_argument(
        "--tp-cols", help="Comma-separated indices for datapoint columns (e.g. 3,4,5)"
    )
    parser.add_argument(
        "--processed", action="store_true", help="Add processed data for the raw table"
    )
    parser.add_argument(
        "--rds", help="Path to the RDS file containing the processed R object"
    )
    args = parser.parse_args()

    if args.csv:
        csv_path = args.csv
    else:
        csv_path = input("Enter path to the CSV file: ").strip()
    if not os.path.exists(csv_path):
        print("File does not exist.")
        return

    open_func = gzip.open if csv_path.endswith(".gz") else open
    with open_func(csv_path, "rt") as f:
        first_line = f.readline()
        if "\t" in first_line:
            sep = "\t"
        else:
            sep = ","
    df = pd.read_csv(
        csv_path, sep=sep, compression="gzip" if csv_path.endswith(".gz") else None
    )
    print(f"CSV loaded with {df.shape[0]} rows and {df.shape[1]} columns.")

    if args.id_cols:
        id_cols = []
        for part in args.id_cols.split(","):
            part = part.strip()
            if ":" in part:
                start, end = part.split(":")
                if start.strip().isdigit() and end.strip().isdigit():
                    id_cols.extend(
                        df.columns[int(i) - 1] for i in range(int(start), int(end) + 1)
                    )
            elif part.isdigit():
                id_cols.append(df.columns[int(part) - 1])
    else:
        id_cols = prompt_column_selection(df, "identifier columns")

    if args.tp_cols:
        tp_cols = []
        for part in args.tp_cols.split(","):
            part = part.strip()
            if ":" in part:
                start, end = part.split(":")
                if start.strip().isdigit() and end.strip().isdigit():
                    tp_cols.extend(
                        df.columns[int(i) - 1] for i in range(int(start), int(end) + 1)
                    )
            elif part.isdigit():
                tp_cols.append(df.columns[int(part) - 1])
    else:
        tp_cols = prompt_column_selection(df, "datapoint columns")

    if args.table:
        table_name = args.table.strip()
    else:
        print(
            "Remember to follow the naming rules: <specie>_<datatype>_<free identifier>_<replicate>."
            "For example 'zebrafish_proteomics_test_1'."
            "For further rules, refer to the README."
        )
        table_name = input("\nEnter a name for the new table in the database: ").strip()
    if not table_name.isidentifier():
        print("Invalid table name.")
        return

    if args.db:
        db_path = args.db
    else:
        db_path = input("\nEnter path to the database: ").strip()
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        df.to_sql(table_name, conn, if_exists="replace", index=False)
        print(f"Table '{table_name}' uploaded to the database.")
    except Exception as e:
        print(f"Error uploading the table to the database: {e}")
        conn.close()
        return

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS table_metadata (
            table_name TEXT PRIMARY KEY,
            identifier_columns TEXT,
            datapoint_columns TEXT
        )
    """)

    id_str = ",".join(id_cols)
    tp_str = ",".join(tp_cols)

    cursor.execute(
        """
        INSERT OR REPLACE INTO table_metadata (table_name, identifier_columns, datapoint_columns)
        VALUES (?, ?, ?)
    """,
        (table_name, id_str, tp_str),
    )
    conn.commit()
    print("Metadata added to 'table_metadata'.")

    # === Step 8: Prompt addition of processed data ===
    add_processed = args.processed

    if add_processed:
        key = f'{table_name}_{"_".join(tp_cols)}_dep'
        if args.rds:
            rds_path = args.rds.strip()
        else:
            rds_path = input(
                "Enter the path to the RDS file containing the processed R object: "
            ).strip()
        if not os.path.exists(rds_path):
            print("RDS file does not exist.")
        else:
            import rpy2.robjects as ro

            ro.r(f"""
            library(storr)
            # Load processed object from RDS
            obj <- readRDS("{rds_path}")
            con <- DBI::dbConnect(RSQLite::SQLite(), "{db_path}")
            # Use connection 'conn' to the database and create cache object
            cache <- storr::storr_dbi(
                tbl_data = "storr_data",
                tbl_keys = "storr_keys",
                con = con
            )      
            # Store object with the generated key
            cache$set(key = "{key}", value = obj)
            """)
            print("Processed data added correctly to cache")
    conn.close()
    print("done")


if __name__ == "__main__":
    main()