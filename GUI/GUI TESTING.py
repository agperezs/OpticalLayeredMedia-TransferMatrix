from tkinter import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

import TMM_Python as TMM
import numpy as np

# Prueba
def graph(wl1,wl2,delta):
    wavelengths = np.arange(float(wl1), float(wl2)+float(delta), float(delta)).tolist()
    fig = Figure(figsize = (10, 2), dpi = 100)
    plot1 = fig.add_subplot(111)
    plot1.plot(wavelengths,wavelengths)
    canvas = FigureCanvasTkAgg(fig, master = root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=1,column=1, columnspan=12)

# Crear la ventana
root = Tk()
root.title("Transfer Matrix Method")
root.geometry("1280x720")
root.iconbitmap("C:/Users/Panda/Downloads/TEC-Copy.ico")

# Agregar un text field
n0_text = Label(root, text= "Ingresar índice de refracción del ambiente:")
n0_text.grid(row=0,column=0,columnspan=2,pady=10)
ns_text = Label(root, text= "Ingresar índice de refracción del sustrato:")
ns_text.grid(row=0,column=3,columnspan=2,pady=10)
wl_text = Label(root, text= "Ingresar rango de longitudes de onda y delta:")
wl_text.grid(row=0,column=6,columnspan=2,padx=1,pady=10)
sep = Label(root, text="-").grid(row=0, column=9)
delta_text = Label(root, text="Delta: ").grid(row=0,column=11)

# Agregar un input field
n0 = Entry(root, width=10, borderwidth=2)
n0.grid(row=0,column=2, padx=10, pady=10)
ns = Entry(root, width=10, borderwidth=2)
ns.grid(row=0,column=5, padx=10, pady=10)
wl1 = Entry(root, width=10, borderwidth=2)
wl1.grid(row=0, column=8)
wl2 = Entry(root, width=10, borderwidth=2)
wl2.grid(row=0, column=10)
delta = Entry(root, width=10, borderwidth=2)
delta.grid(row=0, column=12)


# Agregar un boton
CalcM = Button(root, text="Calcular Matriz de Transferencia", padx=50,pady=20,command=lambda: graph(wl1.get(),wl2.get(),delta.get()))
CalcM.grid(row=1,column=0)
# Crear el loop para la actualizacion del GUI
root.mainloop() 