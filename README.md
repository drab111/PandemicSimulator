# Pandemic Simulator

A simplified pandemic simulation model that visualizes the spread of an infectious disease across four global regions: **Eurasia**, **Africa**, **America**, and **Oceania**.

The simulation is based on an extended **SIRD** model (Susceptible, Infected, Recovered, Dead), with additional support for:
- **Migration** between regions,
- **Demographic changes** (births and unrelated deaths),
- **Real-time visualization** of infected population on a world map.

---

## Features

- Define initial parameters like infection rate, recovery rate, and death rate for each region.
- Start the outbreak in a selected region.
- Run the simulation for a custom number of days.
- Visualize infected population as live color changes on a simplified world map.
- See how different parameters affect the global outcome.

---

## Requirements

Make sure you have the following installed:

- **Python 3.8+**
- **MATLAB** (with Python engine API)
- Python libraries:
  ```bash
  pip install matplotlib numpy cartopy
  ```

> 💡 MATLAB is used to solve differential equations numerically. Python is responsible for visualization.

---

## How to Run

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/pandemic-simulator.git
   cd pandemic-simulator
   ```

2. Open the main script in Python:
   ```
   PandemicSimulator.py
   ```

3. Run the script:
   ```bash
   python PandemicSimulator.py
   ```

4. You will be asked to enter simulation parameters:
   - Starting region
   - Initial number of infected
   - Simulation duration (days)
   - Infection, recovery, and death rates per region
   - Birth and background death rates

5. Watch the simulation live:
   - Infection levels are shown with color gradients on a world map.
   - Time-series plots for each region are also displayed.

---

## File Overview

- `PandemicSimulator.py` – Python interface and visualization
- `symulacja_4kont.m` – MATLAB simulation core

---

## Technologies Used

- **Python** (visualization, user interface)
- **MATLAB** (numerical simulation – SIRD with migration)
- **Cartopy** (geographical mapping)
- **Matplotlib** (plotting)

---

## Notes

- This is a basic educational tool and does **not** represent real pandemic behavior.
- Results are purely illustrative and meant for experimentation and learning.

---

## Contact

Created by **Wiktor Drab**  
GitHub: [drab111](https://github.com/drab111)
