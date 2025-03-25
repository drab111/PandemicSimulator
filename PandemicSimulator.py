import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

import matlab.engine
import numpy as np

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from shapely.ops import unary_union
from matplotlib.patches import Polygon


class PandemicGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Symulator Pandemii")

        # Pola do wprowadzania parametrów:
        # startCont (1,2,3,4)
        ttk.Label(master, text="Kontynent startowy:").grid(row=0, column=0, padx=5, pady=5)
        self.continent_names = ["Eurazja (EA)", "Afryka (AF)", "Ameryka (AM)", "Oceania (OC)"]
        self.combobox_start = ttk.Combobox(master, values=self.continent_names, state="readonly", width=19)
        self.combobox_start.set("Eurazja (EA)")  # Domyślna wartość
        self.combobox_start.grid(row=0, column=1)

        # Mapa nazw kontynentów na liczby
        self.continent_map = {
            "Eurazja (EA)": 1,
            "Afryka (AF)": 2,
            "Ameryka (AM)": 3,
            "Oceania (OC)": 4
        }

        # Początkowa liczba chorych (tylko int)
        ttk.Label(master, text="Początkowa liczba chorych:").grid(row=1, column=0, padx=5, pady=5)
        validate_cmd = (master.register(self.validate_int), '%P')
        self.entry_initial_infected = ttk.Entry(master, validate='key', validatecommand=validate_cmd)
        self.entry_initial_infected.insert(0, "10000")
        self.entry_initial_infected.grid(row=1, column=1)

        # Dni symulacji
        ttk.Label(master, text="Czas symulacji (dni):").grid(row=2, column=0, padx=5, pady=5)
        self.entry_tEnd = ttk.Entry(master)
        self.entry_tEnd.insert(0, "365")
        self.entry_tEnd.grid(row=2, column=1)

        # Beta
        ttk.Label(master, text="Szybkość infekcji (EA AF AM OC):").grid(row=3, column=0, padx=5, pady=5)
        self.entry_beta = ttk.Entry(master)
        self.entry_beta.insert(0, "0.25 0.25 0.25 0.25")
        self.entry_beta.grid(row=3, column=1)

        # Gamma
        ttk.Label(master, text="Szybkość wyzdrowień (EA AF AM OC):").grid(row=4, column=0, padx=5, pady=5)
        self.entry_gamma = ttk.Entry(master)
        self.entry_gamma.insert(0, "0.05 0.05 0.05 0.05")
        self.entry_gamma.grid(row=4, column=1)

        # Mu
        ttk.Label(master, text="Śmiertelność (EA AF AM OC):").grid(row=5, column=0, padx=5, pady=5)
        self.entry_mu = ttk.Entry(master)
        self.entry_mu.insert(0, "0.01 0.01 0.01 0.01")
        self.entry_mu.grid(row=5, column=1)

        # Szybkość narodzin
        ttk.Label(master, text="Szybkość narodzin (EA AF AM OC):").grid(row=6, column=0, padx=5, pady=5)
        self.entry_birth = ttk.Entry(master)
        self.entry_birth.insert(0, "0.0001 0.0001 0.0001 0.0001")
        self.entry_birth.grid(row=6, column=1)

        # Śmiertelność ogólna (nie-pandemiczna)
        ttk.Label(master, text="Ogólna śmiertelność (EA AF AM OC):").grid(row=7, column=0, padx=5, pady=5)
        self.entry_other_death = ttk.Entry(master)
        self.entry_other_death.insert(0, "0.0001 0.0001 0.0001 0.0001")
        self.entry_other_death.grid(row=7, column=1)

        # Przycisk startu symulacji
        self.run_button = ttk.Button(master, text="Uruchom symulację", command=self.run_simulation)
        self.run_button.grid(row=8, column=0, columnspan=2, pady=10)

        # Silnik MATLAB - startowany przy pierwszym użyciu
        self.eng = None

        # Przygotowujemy "geometrie" kontynentów (4 kategorie)
        self.continent_geoms = self.load_continent_polygons()
        # Nazwy w takiej kolejności jak w skrypcie symulacja_4kont (kolumny)
        self.continent_keys = ["Eurazja", "Afryka", "Ameryka", "Oceania"]



    # Walidacja danych wprowadzonych przez użytkownika:
    def validate_int(self, new_value):
        return new_value.isdigit()  # True jeśli same cyfry lub puste
    
    def validate_inputs(self):
        # 1) Liczba chorych > 0
        try:
            infected_val = int(self.entry_initial_infected.get())
            if infected_val < 0:
                raise ValueError("Ujemna liczba chorych")
        except ValueError:
            messagebox.showerror("Błąd danych", "Początkowa liczba chorych musi być nieujemną liczbą całkowitą.")
            return False

        # 2) Ilość symulowanych dni > 0
        try:
            t_end_val = float(self.entry_tEnd.get())
            if t_end_val <= 0:
                raise ValueError("Czas symulacji <= 0")
        except ValueError:
            messagebox.showerror("Błąd danych", "Czas symulacji musi być liczbą dodatnią.")
            return False

        # Sprawdzamy wszystkie współczynniki (muszą po 4 floaty z przedziału [0..1])
        if not self.validate_four_floats(self.entry_beta.get()):
            messagebox.showerror("Błąd danych", "Szybkość infekcji musi zawierać dokładnie 4 wartości typu float z przedziału [0..1].")
            return False

        if not self.validate_four_floats(self.entry_gamma.get()):
            messagebox.showerror("Błąd danych", "Szybkość wyzdrowień musi zawierać dokładnie 4 wartości typu float z przedziału [0..1].")
            return False

        if not self.validate_four_floats(self.entry_mu.get()):
            messagebox.showerror("Błąd danych", "Śmiertelność (pandemiczna) musi zawierać dokładnie 4 wartości typu float z przedziału [0..1].")
            return False
        
        if not self.validate_four_floats(self.entry_birth.get()):
            messagebox.showerror("Błąd danych", "Szybkość narodzin musi zawierać dokładnie 4 wartości typu float z przedziału [0..1].")
            return False

        if not self.validate_four_floats(self.entry_other_death.get()):
            messagebox.showerror("Błąd danych", "Ogólna śmiertelność (inne zgony) musi zawierać dokładnie 4 wartości typu float z przedziału [0..1].")
            return False
        
        return True # gdy wszystko OK

    def validate_four_floats(self, text_value):
        parts = text_value.split()
        if len(parts) != 4:
            return False
        for p in parts:
            try:
                value = float(p)
                if not (0 <= value <= 1):
                    return False
            except ValueError:
                return False
        return True



    # Ładowanie kształtów świata z Natural Earth i podział na 4 kontynenty
    def load_continent_polygons(self):
        shpfilename = shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')

        continents = {
            "Eurazja": [],
            "Afryka": [],
            "Ameryka": [],
            "Oceania": []
        }

        def classify_continent(cont_str):
            if cont_str in ["Europe", "Asia"]:
                return "Eurazja"
            elif cont_str == "Africa":
                return "Afryka"
            elif cont_str in ["North America", "South America"]:
                return "Ameryka"
            elif cont_str == "Oceania":
                return "Oceania"
            return None

        # Mapujemy kontynenty krajów wedle naszego podziału na 4 kontynenty
        for record in shpreader.Reader(shpfilename).records():
            c = record.attributes['CONTINENT']
            group = classify_continent(c)
            if group is not None:
                continents[group].append(record.geometry)

        for k in continents:
            if len(continents[k]) > 0:
                # Sklejamy poszczególne figury reprezentujące kraje w jeden kontynent
                continents[k] = unary_union(continents[k])
            else:
                continents[k] = None

        return continents



    # Funkcja wywoływana po kliknięciu "Uruchom symulację"
    def run_simulation(self):
        # Najpierw sprawdzamy poprawność danych
        if not self.validate_inputs():
            return
        
        # Start MATLAB Engine (jeśli nie wystartowany)
        if not self.eng:
            self.eng = matlab.engine.start_matlab()

        # Czytamy parametry wprowadzone przez użytkownika
        startCont = int(self.continent_map[self.combobox_start.get()])
        tEnd = float(self.entry_tEnd.get())

        beta_list = list(map(float, self.entry_beta.get().split()))
        gamma_list = list(map(float, self.entry_gamma.get().split()))
        mu_list = list(map(float, self.entry_mu.get().split()))
        birth_list = list(map(float, self.entry_birth.get().split()))
        other_death_list = list(map(float, self.entry_other_death.get().split()))
        
        initial_infected = int(self.entry_initial_infected.get())

        # Budujemy strukturę param w formacie zgodnym z MATLAB
        param = {}
        param['startCont'] = float(startCont)
        param['tEnd'] = float(tEnd)
        param['beta'] = matlab.double(beta_list)
        param['gamma'] = matlab.double(gamma_list)
        param['mu'] = matlab.double(mu_list)
        param['birth'] = matlab.double(birth_list)
        param['otherDeath'] = matlab.double(other_death_list)
        param['userInf'] = float(initial_infected)

        # Wywołujemy symulacja_4kont.m
        t, S, I, R, D = self.eng.symulacja_4kont(param, nargout = 5)

        # Konwertujemy wyniki z MATLAB do numpy
        t_py = np.array(t)
        S_py = np.array(S)
        I_py = np.array(I)
        R_py = np.array(R)
        D_py = np.array(D)

        # Animacja mapy
        self.show_map_animation_cartopy(t_py, S_py, I_py, R_py, D_py)


    # Animacja cartopy
    def show_map_animation_cartopy(self, t_py, S_py, I_py, R_py, D_py):
        import cartopy.feature as cfeature

        fig = plt.figure(figsize=(9,6))
        ax_map = plt.axes(projection=ccrs.Robinson())
        ax_map.set_global()
        ax_map.coastlines('110m')
        ax_map.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax_map.set_title("Symulacja pandemii")

        # Mapa kontynentów -> patches
        patches = {}
        for ckey in self.continent_keys:
            geom = self.continent_geoms[ckey]
            if geom is None:
                continue

            if geom.geom_type == 'Polygon':
                geoms = [geom]
            else:
                geoms = list(geom.geoms)

            poly_list = []
            for g in geoms:
                lon, lat = g.exterior.coords.xy
                poly = Polygon(
                    np.array([lon, lat]).T, 
                    facecolor='white', 
                    edgecolor='none', 
                    transform=ccrs.PlateCarree(),
                    alpha=0.7
                )
                ax_map.add_patch(poly)
                poly_list.append(poly)

            patches[ckey] = poly_list

        # Dodajemy oś tekstową na dole figury aby wyświetlać informacje
        ax_text = fig.add_axes([0.05, 0.0, 0.9, 0.15])
        ax_text.axis('off')
        
        # Przygotowujemy teksty dla każdego z 4 kontynentów (4 linie)
        text_lines = []
        for i, ckey in enumerate(self.continent_keys):
            y_pos = 0.8 - 0.25*i
            tline = ax_text.text(
                0.01, y_pos, 
                f"{ckey}: ???", 
                transform=ax_text.transAxes, 
                fontweight='bold', 
                fontsize=10
            )
            text_lines.append(tline)

        plt.ion() 
        plt.show() 

        num_steps = len(t_py)

        # Pętla animacji
        for step in range(num_steps):
            S_now = S_py[step,:]
            I_now = I_py[step,:]
            R_now = R_py[step,:]
            D_now = D_py[step,:]

            N_now = S_now + I_now + R_now
            infected_ratio = np.where(N_now>0, I_now/N_now, 0.0)

            for i, ckey in enumerate(self.continent_keys):
                if ckey not in patches:
                    continue

                color_val = max(0.0, min(1.0, infected_ratio[i]))
                r = color_val
                g = 1 - color_val
                b = 0

                for poly in patches[ckey]:
                    poly.set_facecolor((r,g,b))

                percent_inf = color_val * 100
                dead_num = D_now[i]
                recov_num = R_now[i]
                text_lines[i].set_text(
                    f"{ckey}: {percent_inf:.2f}% zarażonych,   Zmarli={dead_num:,.0f},   Wyleczeni={recov_num:,.0f}"
                )

            plt.pause(0.25)

        plt.ioff()
        plt.show()

        last_step = num_steps - 1
        S_end = S_py[last_step,:]
        I_end = I_py[last_step,:]
        R_end = R_py[last_step,:]
        D_end = D_py[last_step,:]

        total_infected = np.sum(I_end)
        total_dead = np.sum(D_end)
        total_recovered = np.sum(R_end)

        summary_msg = (
            f"Łączna liczba osób zmarłych: {total_dead:,.0f}\n"
            f"Łączna liczba osób wyleczonych: {total_recovered:,.0f}\n"
            f"Łączna liczba osób aktualnie zainfekowanych: {total_infected:,.0f}"
        )
        messagebox.showinfo("Wyniki końcowe symulacji", summary_msg)



# Main
if __name__ == "__main__":
    root = tk.Tk()
    app = PandemicGUI(root)
    root.mainloop()