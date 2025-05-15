# FastaComparer

## Generator Sekwencji FASTA

Prosty program w Pythonie do generowania losowych sekwencji DNA w formacie FASTA. Umożliwia personalizację (ID, opis, wstawienie imienia), oblicza podstawowe statystyki sekwencji oraz oferuje opcjonalne porównanie wygenerowanej sekwencji z fragmentem sekwencji referencyjnej pobranej z bazy NCBI.

### Funkcjonalności

*   Generowanie losowych sekwencji DNA (A, C, G, T) o zadanej przez użytkownika długości.
*   Możliwość zdefiniowania ID i opisu sekwencji, które są umieszczane w nagłówku pliku FASTA.
*   Opcjonalne wstawienie podanego imienia w losowe miejsce wygenerowanej sekwencji (imię nie wpływa na statystyki ani długość liczoną do statystyk).
*   Zapis wygenerowanej sekwencji do pliku `.fasta` z nazwą odpowiadającą podanemu ID.
*   Formatowanie sekwencji w pliku FASTA z zawijaniem linii dla lepszej czytelności i zgodności ze standardem.
*   Wyświetlanie statystyk wygenerowanej sekwencji:
    *   Procentowa zawartość każdego nukleotydu (A, C, G, T).
    *   Procentowa zawartość par C-G.
*   **Nowość:** Możliwość pobrania sekwencji referencyjnej z NCBI (np. chromosomu ludzkiego) za pomocą jej ID dostępu.
*   **Nowość:** Porównanie wygenerowanej losowo sekwencji z początkowym fragmentem pobranej sekwencji referencyjnej i wyświetlenie procentu zgodności.

### Wymagania

*   Python 3.x
*   Biblioteka `requests` (do pobierania danych z NCBI)

### Instalacja biblioteki `requests`

Jeśli nie masz zainstalowanej biblioteki `requests`, możesz ją zainstalować za pomocą pip:

```bash
pip install requests
