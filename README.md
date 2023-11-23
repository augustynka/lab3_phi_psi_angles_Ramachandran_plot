# lab3_zad2b



Program wykorzystuje biblioteki BioPython (PDBList, PDBParser, DSSP, PPBuilder), matplotlib.pyplot, seaborn, pandas, numpy

_Funkcje:_

**calculate_phi_psi_angles(structure, dssp_path):**
na podstawie przekazanej struktury białka i pliku dssp funkcja buduje peptydy dla każdego łańcucha w modelu, następnie pobiera listę kątów phi psi, jeśli istnieją oba, pobiera strukturę drugorzędową z DSSP i dodaje dane (phi, psi, struktura drugorzędowa) do tablicy danych phi_psi_data.

**plot_ramachandran(phi_psi_data):**

tworzy ramkę danych przypisującą każdej reszcie aminokwasowej kąty i strukturę drugorzędową

tworzy scatterplot, na którym kolory kropek odpowiadają wskazanym strukturom drugorzędowym.

zapisuje wykres do pliku .png

_Działanie programu:_

dla wskazanego w kodzie pdb_id program pobiera plik z bazy PDB, następnie załadowuje z niego strukturę. 

Obliczane są kąty phi, psi za pomocą funkcji calculate_phi_psi_angles(structure,’sciezka_do_dssp’)

Tworzony jest wykres Ramachandrana zapisany do pliku .png.
