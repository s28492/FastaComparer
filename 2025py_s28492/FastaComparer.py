# Cel programu: Generowanie losowych sekwencji DNA w formacie FASTA
#               z możliwością personalizacji (ID, opis, wstawione imię),
#               wyświetlanie podstawowych statystyk sekwencji,
#               oraz opcjonalne porównanie z sekwencją referencyjną z NCBI.
# Kontekst zastosowania: Proste narzędzie edukacyjne, do generowania
#                        danych testowych lub demonstracji użycia API bioinformatycznych.

import random  # Import modułu random do generowania losowych wartości
import re  # Import modułu re do obsługi wyrażeń regularnych (dla walidacji ID)
import requests  # Import modułu requests do komunikacji z API NCBI


# ORIGINAL:
# def generuj_sekwencje_dna(dlugosc):
#     # Funkcja generuje losową sekwencję DNA o podanej długości.
#     # Argumenty:
#     #   dlugosc (int): Pożądana długość sekwencji DNA.
#     # Zwraca:
#     #   str: Losowo wygenerowana sekwencja DNA.
#     nukleotydy = ['A', 'C', 'G', 'T'] # Lista możliwych nukleotydów
#     sekwencja = ''.join(random.choice(nukleotydy) for _ in range(dlugosc)) # Generowanie sekwencji
#     return sekwencja
# MODIFIED (bez zmian w logice, tylko dodano komentarz dla jasności):
def generuj_sekwencje_dna(dlugosc):
    # Funkcja generuje losową sekwencję DNA o podanej długości.
    # Argumenty:
    #   dlugosc (int): Pożądana długość sekwencji DNA.
    # Zwraca:
    #   str: Losowo wygenerowana sekwencja DNA.
    nukleotydy = ['A', 'C', 'G', 'T']  # Lista możliwych nukleotydów
    sekwencja = ''.join(random.choice(nukleotydy) for _ in range(dlugosc))  # Generowanie sekwencji
    return sekwencja


def oblicz_statystyki(sekwencja_dna_oryginalna):
    # Funkcja oblicza statystyki nukleotydowe dla podanej sekwencji DNA.
    # Ważne: ta funkcja powinna otrzymać CZYSTĄ sekwencję DNA, bez wstawionego imienia.
    # Argumenty:
    #   sekwencja_dna_oryginalna (str): Sekwencja DNA, dla której liczone są statystyki.
    # Zwraca:
    #   dict: Słownik zawierający procentową zawartość A, C, G, T oraz %CG.
    dlugosc_oryginalna = len(sekwencja_dna_oryginalna)  # Długość oryginalnej sekwencji
    if dlugosc_oryginalna == 0:  # Zabezpieczenie przed dzieleniem przez zero
        return {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'CG': 0.0}

    licznik_A = sekwencja_dna_oryginalna.count('A')  # Liczba nukleotydów A
    licznik_C = sekwencja_dna_oryginalna.count('C')  # Liczba nukleotydów C
    licznik_G = sekwencja_dna_oryginalna.count('G')  # Liczba nukleotydów G
    licznik_T = sekwencja_dna_oryginalna.count('T')  # Liczba nukleotydów T

    # Obliczanie procentowej zawartości
    procent_A = (licznik_A / dlugosc_oryginalna) * 100
    procent_C = (licznik_C / dlugosc_oryginalna) * 100
    procent_G = (licznik_G / dlugosc_oryginalna) * 100
    procent_T = (licznik_T / dlugosc_oryginalna) * 100
    procent_CG = ((licznik_C + licznik_G) / dlugosc_oryginalna) * 100

    return {  # Zwracanie słownika ze statystykami
        'A': procent_A,
        'C': procent_C,
        'G': procent_G,
        'T': procent_T,
        'CG': procent_CG
    }


def wstaw_imie_do_sekwencji(sekwencja, imie):
    # Funkcja wstawia podane imię w losowe miejsce w sekwencji.
    # Argumenty:
    #   sekwencja (str): Sekwencja DNA, do której ma być wstawione imię.
    #   imie (str): Imię do wstawienia.
    # Zwraca:
    #   str: Sekwencja DNA z wstawionym imieniem.
    if not imie:  # Jeśli imię jest puste, nie wstawiaj niczego
        return sekwencja
    if not sekwencja:  # Jeśli sekwencja jest pusta, imię jest wstawiane na początku
        return imie
    # ORIGINAL:
    # pozycja_wstawienia = random.randint(0, len(sekwencja)) # Losowanie pozycji wstawienia
    # MODIFIED (zmiana losowania pozycji, aby imię nie było wstawiane na samym końcu, co może wyglądać nienaturalnie bez następujących nukleotydów):
    # Uzasadnienie: Poprawia estetykę wstawienia, unikając sytuacji, gdzie imię jest ostatnim elementem sekwencji bez nukleotydów po nim,
    # chyba że oryginalna sekwencja była pusta lub bardzo krótka.
    if len(sekwencja) > 0:
        pozycja_wstawienia = random.randint(0, len(sekwencja) - 1 if len(sekwencja) > 0 else 0)
    else:
        pozycja_wstawienia = 0

    sekwencja_z_imieniem = sekwencja[:pozycja_wstawienia] + imie + sekwencja[
                                                                   pozycja_wstawienia:]  # Tworzenie nowej sekwencji
    return sekwencja_z_imieniem


# NOWA FUNKCJA (Ulepszenie 3: Formatowanie FASTA z zawijaniem linii)
def formatuj_sekwencje_do_zapisu_fasta(sekwencja, dlugosc_linii=70):
    # Funkcja formatuje długą sekwencję na wiele linii o określonej długości.
    # Argumenty:
    #   sekwencja (str): Sekwencja do sformatowania.
    #   dlugosc_linii (int): Maksymalna długość każdej linii sekwencji.
    # Zwraca:
    #   str: Sekwencja sformatowana z przejściami do nowej linii.
    # Uzasadnienie: Standard FASTA często używa zawijania linii dla czytelności i kompatybilności.
    if not sekwencja:
        return ""
    return "\n".join(sekwencja[i:i + dlugosc_linii] for i in range(0, len(sekwencja), dlugosc_linii))


# NOWA FUNKCJA (Ulepszenie 1: Pobieranie sekwencji z NCBI)
def pobierz_sekwencje_z_ncbi(accession_id):
    # Funkcja pobiera sekwencję w formacie FASTA z NCBI E-utilities.
    # Argumenty:
    #   accession_id (str): Identyfikator dostępu NCBI (np. "NC_000024.10").
    # Zwraca:
    #   str: Czysta sekwencja DNA (bez nagłówka FASTA i znaków nowej linii) lub None w przypadku błędu.
    # Uzasadnienie: Umożliwia porównanie generowanej sekwencji z rzeczywistymi danymi genomowymi.
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession_id,
        "rettype": "fasta",
        "retmode": "text"
    }
    print(f"\nPobieranie sekwencji {accession_id} z NCBI...")
    try:
        response = requests.get(base_url, params=params, timeout=30)  # Timeout w przypadku za długiego działania
        response.raise_for_status()  # Wyjątek dla błędów HTTP

        fasta_content = response.text
        # Usunięcie nagłówka FASTA i znaków nowej linii, aby uzyskać czystą sekwencję
        lines = fasta_content.splitlines()
        czysta_sekwencja = "".join(lines[1:])  # Pomijamy nagłóek
        # Usunięcie innych znaków
        czysta_sekwencja = re.sub(r'[^ATGCatgc]', '', czysta_sekwencja).upper()
        print(f"Pobrano {len(czysta_sekwencja)} nukleotydów.")
        return czysta_sekwencja
    except requests.exceptions.HTTPError as http_err:
        print(f"Błąd HTTP podczas pobierania z NCBI: {http_err}")
        print(f"Treść odpowiedzi serwera (fragment): {response.text[:500]}...")
    except requests.exceptions.RequestException as req_err:
        print(f"Błąd żądania podczas pobierania z NCBI: {req_err}")
    except Exception as e:
        print(f"Nieoczekiwany błąd podczas pobierania z NCBI: {e}")
    return None


# NOWA FUNKCJA (Ulepszenie 1: Porównywanie sekwencji)
def porownaj_sekwencje(gen_seq, ref_seq_fragment):
    # Funkcja porównuje dwie sekwencje nukleotyd po nukleotydzie.
    # Zakłada, że obie sekwencje mają tę samą długość.
    # Argumenty:
    #   gen_seq (str): Wygenerowana sekwencja.
    #   ref_seq_fragment (str): Fragment sekwencji referencyjnej.
    # Zwraca:
    #   float: Procent zgodności.
    # Uzasadnienie: Pozwala ocenić, na ile losowa sekwencja odpowiada fragmentowi rzeczywistej.
    matches = 0
    dlugosc = len(gen_seq)  # Obie powinny mieć tę samą długość do porównania
    if dlugosc == 0:
        return 0.0
    for i in range(dlugosc):
        if gen_seq[i] == ref_seq_fragment[i]:
            matches += 1
    return (matches / dlugosc) * 100


def main():
    # Wymaganie 2: Długość sekwencji określana przez użytkownika
    while True:
        try:
            dlugosc_sek = int(input("Podaj długość sekwencji (np. 20-1000): "))
            if dlugosc_sek <= 0:  # Zmieniono z < 0 na <=0, bo długość 0 nie ma sensu
                print("Długość sekwencji musi być liczbą dodatnią. Spróbuj ponownie.")
                continue
            break
        except ValueError:
            print("Nieprawidłowa wartość. Długość musi być liczbą całkowitą.")

    # Wymaganie 3: Zapytanie o nazwę (ID) i opis sekwencji
    # ORIGINAL (dla ID):
    # id_sek = input("Podaj ID sekwencji: ")
    # MODIFIED (Ulepszenie 2: Walidacja ID sekwencji):
    # Uzasadnienie: Zapobiega błędom przy tworzeniu plików z nieprawidłowymi nazwami.
    while True:
        id_sek = input("Podaj ID sekwencji (tylko litery, cyfry, '_' lub '-'): ")
        if re.match(r"^[a-zA-Z0-9_-]+$", id_sek):  # Sprawdzenie poprawności ID
            break
        else:
            print("Nieprawidłowy ID. Użyj tylko liter, cyfr, podkreślenia (_) lub myślnika (-).")

    opis_sek = input("Podaj opis sekwencji: ")

    # Wymaganie 6 (część): Zapytanie o imię do wstawienia
    imie_do_wstawienia = input("Podaj imię (pozostaw puste, jeśli nie chcesz wstawiać): ")

    # Wymaganie 1: Generowanie losowej sekwencji DNA
    oryginalna_sekwencja_dna = generuj_sekwencje_dna(dlugosc_sek)

    # Wymaganie 6: Wstawienie imienia w losowe miejsce
    sekwencja_z_imieniem_do_pliku = wstaw_imie_do_sekwencji(oryginalna_sekwencja_dna, imie_do_wstawienia)

    # Wymaganie 4: Zapis wyniku do pliku FASTA
    nazwa_pliku = f"{id_sek}.fasta"
    naglowek_fasta = f">{id_sek} {opis_sek}"

    try:
        with open(nazwa_pliku, "w") as plik:
            plik.write(naglowek_fasta + "\n")
            # ORIGINAL (zapis sekwencji):
            # plik.write(sekwencja_z_imieniem_do_pliku + "\n")
            # MODIFIED (Ulepszenie 3: Formatowanie FASTA z zawijaniem linii):
            # Uzasadnienie: Zgodność ze standardem FASTA i lepsza czytelność.
            sformatowana_sekwencja_do_zapisu = formatuj_sekwencje_do_zapisu_fasta(sekwencja_z_imieniem_do_pliku, 70)
            plik.write(sformatowana_sekwencja_do_zapisu + "\n")
        print(f"\nSekwencja została zapisana do pliku {nazwa_pliku}")
    except IOError as e:
        print(f"Błąd zapisu do pliku {nazwa_pliku}: {e}")
        return

    # Wymaganie 5: Wyświetlanie statystyk sekwencji (na podstawie oryginalnej sekwencji)
    statystyki = oblicz_statystyki(oryginalna_sekwencja_dna)

    print("\nStatystyki wygenerowanej sekwencji (na podstawie oryginalnej, bez imienia):")
    print(f"A: {statystyki['A']:.1f}%")
    print(f"C: {statystyki['C']:.1f}%")
    print(f"G: {statystyki['G']:.1f}%")
    print(f"T: {statystyki['T']:.1f}%")
    print(f"%CG: {statystyki['CG']:.1f}%")

    # Wyświetlenie przykładowej zawartości pliku
    print(f"\nPrzykładowa zawartość pliku {nazwa_pliku}:")
    print('"""')
    print(naglowek_fasta)
    print(formatuj_sekwencje_do_zapisu_fasta(sekwencja_z_imieniem_do_pliku, 70))  # Również sformatowana
    print('"""')

    # Ulepszenie 1: Porównanie z sekwencją z NCBI
    czy_porownac = input(
        "\nCzy chcesz porównać wygenerowaną sekwencję z sekwencją referencyjną z NCBI? (tak/nie): ").strip().lower()
    if czy_porownac == 'tak':
        ncbi_id = input("Podaj ID dostępu sekwencji z NCBI (np. NC_000024.10 dla ludzkiego chromosomu Y): ").strip()
        if ncbi_id:
            sekwencja_referencyjna = pobierz_sekwencje_z_ncbi(ncbi_id)
            if sekwencja_referencyjna:
                if len(sekwencja_referencyjna) < dlugosc_sek:
                    print(
                        f"Pobrana sekwencja referencyjna ({len(sekwencja_referencyjna)} bp) jest krótsza niż wygenerowana ({dlugosc_sek} bp). Porównanie nie jest możliwe w tej formie.")
                else:
                    # Bierzemy fragment sekwencji referencyjnej o tej samej długości co nasza wygenerowana
                    fragment_referencyjny = sekwencja_referencyjna[:dlugosc_sek]

                    # Porównujemy ORYGINALNĄ wygenerowaną sekwencję (bez imienia)
                    procent_zgodnosci = porownaj_sekwencje(oryginalna_sekwencja_dna, fragment_referencyjny)
                    print(
                        f"\nPorównanie wygenerowanej sekwencji (długość: {dlugosc_sek} bp) z początkowym fragmentem sekwencji {ncbi_id}:")
                    print(f"Procent zgodności nukleotydów: {procent_zgodnosci:.2f}%")
            else:
                print("Nie udało się pobrać sekwencji referencyjnej.")
        else:
            print("Nie podano ID sekwencji z NCBI.")


if __name__ == "__main__":
    print("Generator Sekwencji FASTA v2.0")
    print("--------------------------------")
    print("Pamiętaj, że do pobierania danych z NCBI potrzebne jest połączenie z internetem.")
    print("Upewnij się, że masz zainstalowaną bibliotekę 'requests': pip install requests\n")
    main()
    print("\nZakończono działanie programu.")