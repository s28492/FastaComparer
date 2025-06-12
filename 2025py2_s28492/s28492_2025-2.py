import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
import time

def run_pipeline():
    email = input("Email do NCBI: ")
    taxid = input("TaxID organizmu: ")
    min_len = int(input("Minimum length of sequence: "))
    max_len = int(input("Maximum length of sequence: "))
    Entrez.email = email
    print(f"Pobieranie dla TaxID: {taxid} (dł: {min_len}-{max_len})...")
    filtered_data = []
    
    handle_search = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]",
                                   usehistory="y", idtype="acc")
    search_results = Entrez.read(handle_search)
    handle_search.close()

    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    print(f"Znaleziono {count} rekordów.")

    if count == 0: 
        return

    # Batch processing
    batch_s = 100
    for start in range(0, count, batch_s):
        handle_fetch = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                     retstart=start, retmax=batch_s,
                                     webenv=webenv, query_key=query_key)
        for record in SeqIO.parse(handle_fetch, "gb"):
            if min_len <= len(record.seq) <= max_len:
                filtered_data.append({
                    "acc": record.id,
                    "len": len(record.seq),
                    "desc": record.description[:60]
                })
        handle_fetch.close()
        time.sleep(0.3)  # Nie ma za co NCBI
    


    if not filtered_data: 
        print("Brak pobranych rekordów.")
        return

    # Tworzenie df i zapisywanie do .csv
    df = pd.DataFrame(filtered_data).sort_values("len", ascending=False)
    base_file_name = f"taxid_{taxid}_{min_len}-{max_len}"
    df.to_csv(f"{base_file_name}_report.csv", index=False)
    print(f"Zapisano raport: {base_file_name}_report.csv")

    # Wygeneruj i zapisz wykres
    plot_df = df.head(50) if len(df) > 50 else df
    title = f"{'Top 50 z ' if len(df) > 50 else ''}{len(df)} sekw. dla {taxid}"

    plt.figure(figsize=(10, 5))  # Rozmiar wykresu
    plt.plot(plot_df["acc"], plot_df["len"], marker='.', linestyle='-')
    plt.xticks(rotation=75, ha="right", fontsize=7) 
    plt.ylabel("Długość sekwencji")
    plt.xlabel("Numer akcesyjny")  
    plt.title(title, fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{base_file_name}_plot.png")
    plt.show()


if __name__ == "__main__":
    run_pipeline()