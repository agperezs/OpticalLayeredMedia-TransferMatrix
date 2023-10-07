###########################################
#  Transfer Matrix Method Implementation  #
#       Armando Gonzalo Pérez Solis       #
#     Last Updated: 17/09/23 ~ 22:00      #
###########################################

## Importar librerias necesarias
import numpy as np
import matplotlib.pyplot as plt
import math

## Crear un dispositivo
class Layer:
    def __init__(self, refractive_index, thickness):
        self.refractive_index = refractive_index
        self.thickness = thickness

def create_DBR(Layer0, LayerH, LayerL, N, LayerS):
    DBR = [LayerH,LayerL]*N
    Device = [Layer0]
    Device.extend(DBR)
    Device.append(LayerS)
    return Device


def snells_law(n1,n2,inc_angle):
    """
    Calcular el angulo de refracción usando la ley de Snell.
    Python utiliza valores en radianes entonces se hacen las conversiones necesarias.

    Args:
    n1: Refractive index of the first medium.
    n2: Refractive index of the second medium.
    inc_angle: Incident angle in degrees.

    Returns:
    refr_angle: Refracted angle in degrees.
    """

    inc_angle_rad = math.radians(inc_angle)
    refr_angle_rad = math.asin((n1/n2)*math.sin(inc_angle_rad))
    refr_angle = math.degrees(refr_angle_rad)
    return refr_angle

def dynamic_matrix(n, thetad):
    """
    Calculo de la matriz dinamica de la capa 
    Args:
    n: Refractive index of the layer
    thetad: Incident angle in degrees.

    Returns:
    Ds,Dp: Dynamic matrix for both polarizations.
    """
    
    # Calcular valor de angulo en radianes
    theta = math.radians(thetad)
    ## Primero calcular los coeficientes de transmision y reflexion
    # TE - s-polarization
    Ds = np.array([[1, 1],[n*math.cos(theta), -n*math.cos(theta)]], dtype=np.complex128)
    # TM - P-polarization
    Dp = np.array([[math.cos(theta), math.cos(theta)],[n, -n]], dtype=np.complex128)

    return Ds, Dp

def propagation_matrix(n,thetad,wavelength_nm,d_nm):
    """
    Calculo de la matriz de propagacion de la capa
    Args:
    n: Refractive index of the medium.
    thetad: Angle in degrees.
    wavelength: Wavelength
    d: Layer width

    Returns:
    P: Propagation matrix.
    """

    # Pasar angulo a radianes y calcular constantes
    theta = math.radians(thetad)
    wavelength = wavelength_nm/1e9
    k0 = (2*np.pi)/wavelength
    k = k0*n*math.cos(theta)
    d = d_nm/1e9
    P = np.array([[(np.exp(1j*k*d)), 0],[0, (np.exp(-1j*k*d))]],dtype=np.complex128)
     
    return P


def transfer_matrix(Ds,Dp,P):
    """
    Calculo de la transfer matrix 
    Args:
    Ds: s-polarized Dynamic Matrix of each layer of the device
    Dp: p-polarized Dynamic Matrix of each layer of the device
    P: Propagation matrix of each layer of the device

    Returns:
    Ms,Mp: Transfer Matrix for that wave length
    """
### PENDIENTE
def total_transfer_matrix(Device, thetad, wavelengths):
    """
    Calculo de la transfer matrix para el sistema multicapa
    Args:
    layers: List of refractive indices and thickness for each layer.
    thetad: Incidence angle in degrees.
    wavelength: Wavelength range

    Returns:
    Ms,Mt: Total transfer matrix for both polarizations
    """
    # Obtener los indices de refraccion del medio incidente y sustrato
    n0 = Device[0].refractive_index
    ns = Device[-1].refractive_index

    # Desde el inicio se puede obtener la matriz dinamica del ambiente
    D0_s,D0_p = dynamic_matrix(n0, thetad)

    # Inicializar las matrices totales
    Ms_total = []
    Mp_total = []
    Ms,Mp = 0

    # Ciclar por cada longitud de onda
    for wavelength in wavelengths:
        # Iniciar valor de angulo en angulo incidente
        theta = thetad

        # Iniciar matrices vacias
        P_wl = np.empty((len(Device), 2, 2), dtype=complex)
        Ds_wl = np.empty((len(Device), 2, 2), dtype=complex)
        Dp_wl = np.empty((len(Device), 2, 2), dtype=complex)

        # Como no se necesita la matriz de propagacion en el ambiente y sustrato se sustituye por la identidad
        P_wl[0] = P_wl[-1] = np.identity(2)

        # Se hace el calculo de la matriz dinamica del ambiente
        Ds_wl[0] = D0_s
        Dp_wl[0] = D0_p

        # Ciclar por cada capa excluyendo medio incidente y sustrato
        for i in range(1, len(Device) - 1):
            n_last = Device[i-1].refractive_index
            n = Device[i].refractive_index
            d = Device[i].refractive_index

            ## VERIFICAR CUAL THETA VA ACA
            theta = snells_law(n_last,n,theta)

            P = propagation_matrix(n,theta,wavelength,d)
            Ds, Dp = dynamic_matrix(n, theta)

            P_wl[i] = P
            Ds_wl[i] = Ds
            Dp_wl[i] = Dp
        


        return Ms, Mp

## PENDIENTE
def system_response (Ms,Mp,n0,ns):
    """
    Calculo de la reflexión y transmisión del sistema
    Args:
    Ms: Complete transfer Matrix (TE)
    Mp: Complete transfer Matrix (TM)
    n0: Refractive index of incident medium
    ns: Refractive index of substrate

    Returns:
    P: Reflection and Transmission.
    """

    ## Primero calcular los coeficientes de transmision y reflexion
    # TE - s-polarization
    rs = (Ms[1,0])/(Ms[0,0])
    ts = 1/(Ms[0,0])
    # TM - P-polarization
    rp = (Mp[1,0])/(Ms[0,0])
    tp = 1/(Mp[0,0])

    ## Calcular reflexion y transmision
    Rs = (abs(rs))**2 
    Rp = (abs(rp))**2 

    Ts = 0
    Tp = 0

    return Rs, Rp, Ts, Tp

## Definicion de caractersiticas (Entradas)
layer0 = Layer(1, None)
layer1_n = Layer(2.5, 176); layer2_n = Layer(1.55, 264); N = 1
layersub = Layer(1.65, None)
wl1 = 500; wl2 = 500; delta = 0.5 # en nm
incident_angle = 45 # en grados

## Crear el dispositivo y la lista de longitudes de onda
Device = create_DBR(layer0, layer1_n, layer2_n, N, layersub)
wavelength = np.arange(wl1, wl2+delta, delta).tolist()
print("Wavelengths: ", wavelength)
for i in range(len(Device)):
            print("Layer: ", i, " - Refractive Index: ", "{:.2f}".format(Device[i].refractive_index), " ~ Thickness: ", Device[i].thickness, "nm")

## TESTING
print("Refractive Index: ", Device[0].refractive_index)
D0s, D0p = dynamic_matrix(Device[0].refractive_index, incident_angle)
print("Ambient")
print(D0s)
print("Layer 1")
theta1 = snells_law(Device[0].refractive_index,Device[1].refractive_index, incident_angle)
print("Theta1: ", theta1)
D1s, D1p = dynamic_matrix(Device[1].refractive_index, theta1)
print("Dynamic Matrix 1")
print(D1s)
P1 = propagation_matrix(Device[1].refractive_index,theta1,wavelength[0],Device[1].thickness)
print("Propagation Matrix1")
print(P1)
print("Layer 2")
theta2 = snells_law(Device[1].refractive_index,Device[2].refractive_index, theta1)
print("Theta2: ", theta2)
D2s, D2p = dynamic_matrix(Device[2].refractive_index, theta2)
print("Dynamic Matrix 2")
print(D2s)
P2 = propagation_matrix(Device[2].refractive_index,theta2,wavelength[0],Device[2].thickness)
print("Propagation Matrix2")
print(P2)
print("Layer 3")
theta3 = snells_law(Device[2].refractive_index,Device[3].refractive_index, theta2)
print("Theta3: ", theta3)
D3s, D3p = dynamic_matrix(Device[3].refractive_index, theta3)
print("Dynamic Matrix 3")
print(D3s)

print("Calculating TMM \n")
D0sinv=np.linalg.inv(D0s);D1sinv=np.linalg.inv(D1s);D2sinv=np.linalg.inv(D2s);
D0pinv=np.linalg.inv(D0p);D1pinv=np.linalg.inv(D1p);D2pinv=np.linalg.inv(D2p);
Ms = D0sinv@D1s@P1@D1sinv@D2s@P2@D2sinv@D3s
Mp = D0pinv@D1p@P1@D1pinv@D2p@P2@D2pinv@D3p
print("Total Transfer Matrix - s pol")
print(Ms)

print("\n Total Transfer Matrix - p pol")
print(Mp)