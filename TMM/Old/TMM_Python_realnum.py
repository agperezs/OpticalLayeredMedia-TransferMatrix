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


def layer_angle_list(Device, inc_angle):
    """
    Calcular una lista de angulos para el dispositivo ya que no cambia por longitud de onda.
    Python utiliza valores en radianes entonces se hacen las conversiones necesarias.

    Args:
    Device: Contains list of all layers and its properties
    inc_angle: Incident angle in degrees.

    Returns:
    angles: Refracted angles of all device degrees.
    """
    # Inicia lista vacia del tamaño del dispositivo
    angles = list(range(len(Device)))

    # Se introduce el valor de angulo de incidencia
    angles[0] = inc_angle
    theta = inc_angle

    # Cicla a traves del resto de capas para obtener el valor de angulo
    for i in range(1, len(Device)):
            # Calcula el angulo refractado y lo guarda en la lista
            theta = snells_law(Device[i-1].refractive_index,Device[i].refractive_index,theta)
            angles[i] = theta
    
    return angles


def inv_mat(Matrix):
    """
    Calculo de la inversa de una matriz 
    Args:
    Matrix: Any matrix

    Returns:
    InvMat: Inverse of matrix
    """

    InvMat = np.linalg.inv(Matrix)
    return InvMat


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

    # Se inicializan matrices con el primer valor a multiplicar:
    Ms = inv_mat(Ds[0])
    Mp = inv_mat(Dp[0])

    # Se cicla por los valores de las capas internas Π(Dl*Pl*Dl^-1)
    for i in range(1, len(Ds) - 1):
         ## For TE Polarization
         invDs = inv_mat(Ds[i])
         Ms = Ms @ Ds[i] @ P[i] @ invDs

         ## For TM Polarization
         invDp = inv_mat(Dp[i])
         Mp = Mp @ Dp[i] @ P[i] @ invDp

    Ms = Ms @ Ds[-1]
    Mp = Mp @ Dp[-1]
    return Ms, Mp


def total_transfer_matrix(Device, angle_list, wavelengths):
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

    # Se calcula la lista de angulos del dispositivo
    angle_list = layer_angle_list(Device, incident_angle)

    # Desde el inicio se puede obtener la matriz dinamica del ambiente y sustrato
    D0_TE,D0_TM = dynamic_matrix(n0, angle_list[0])
    Ds_TE,Ds_TM = dynamic_matrix(ns, angle_list[-1])

    # Inicializar las matrices totales
    Ms_total = np.empty((len(wavelengths), 2, 2), dtype=complex)
    Mp_total = np.empty((len(wavelengths), 2, 2), dtype=complex)

    # Ciclar por cada longitud de onda
    for e, wavelength in enumerate(wavelengths):

        # Iniciar un array de matrices del tamaño del dispositivo
        P_wl = np.empty((len(Device), 2, 2), dtype=complex)
        Ds_wl = np.empty((len(Device), 2, 2), dtype=complex)
        Dp_wl = np.empty((len(Device), 2, 2), dtype=complex)

        # Como no se necesita la matriz de propagacion en el ambiente y sustrato se sustituye por la identidad
        P_wl[0] = P_wl[-1] = np.identity(2)

        # Se introduce en la lista los valores de matrices dinamicas calculadas anteriormente
        Ds_wl[0] = D0_TE
        Dp_wl[0] = D0_TM
        Ds_wl[-1] = Ds_TE
        Dp_wl[-1] = Ds_TM


        # Ciclar por cada capa excluyendo medio incidente y sustrato
        for i in range(1, len(Device) - 1):
            # Se extrae las caracteristicas de la capa
            n = Device[i].refractive_index
            d = Device[i].thickness
            theta_n = angle_list[i]

            # Se calculan las matrices necesarias
            P = propagation_matrix(n,theta_n,wavelength,d)
            Ds, Dp = dynamic_matrix(n, theta_n)

            # Se guardan en la lista total
            P_wl[i] = P
            Ds_wl[i] = Ds
            Dp_wl[i] = Dp

        Ms, Mp = transfer_matrix(Ds_wl,Dp_wl,P_wl)
        Ms_total[e] = Ms
        Mp_total[e] = Mp
        
    return Ms_total, Mp_total

## PENDIENTE
def system_response(Ms_total,Mp_total,n0,ns,angles):
    """
    Calculo de la reflexión y transmisión del sistema
    Args:
    Ms: Complete transfer Matrix (TE) for all wavelengths
    Mp: Complete transfer Matrix (TM) for all wavelengths
    n0: Refractive index of incident medium
    ns: Refractive index of substrate
    angles: List of angles of system

    Returns:
    P: Reflection and Transmission for al wavelengths.
    """
    # Calcular valores de angulos importantes
    theta0 = math.radians(angles[0])
    thetas = math.radians(angles[-1])

    k_transmittance = (ns*math.cos(thetas))/(n0*math.cos(theta0))

    # Se calcula los valores para todas las longitudes de onda
    for (Ms,Mp) in zip(Ms_total,Mp_total):
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

        Ts = k_transmittance*(abs(ts))**2 
        Tp = k_transmittance*(abs(tp))**2 

    return Rs,Ts,Rp,Tp

## Definicion de caractersiticas (Entradas)
layer0 = Layer(1, None)
layer1_n = Layer(2.5, 176); layer2_n = Layer(1.55, 264); N = 1
layersub = Layer(1.65, None)
wl1 = 500; wl2 = 500; delta = 0.5 # en nm
incident_angle = 45 # en grados

## Crear el dispositivo y la lista de longitudes de onda
Device = create_DBR(layer0, layer1_n, layer2_n, N, layersub)
wavelengths = np.arange(wl1, wl2+delta, delta).tolist()
print("Wavelengths: ", wavelengths)

## Mostrar dispositivo 
#for i in range(len(Device)):
            #print("Layer: ", i, " - Refractive Index: ", "{:.2f}".format(Device[i].refractive_index), " ~ Thickness: ", Device[i].thickness, "nm")

angle_list = layer_angle_list(Device,incident_angle)

Ms_total,Mp_total = total_transfer_matrix(Device, angle_list, wavelengths)
#print("Ms is:")
#print(Ms_total)
#print("\n Mp is:")
#print(Mp_total)

Rs,Ts,Rp,Tp=system_response(Ms_total,Mp_total,Device[0].refractive_index,Device[-1].refractive_index,angle_list)
#print("s - polarized")
print("System Reflectance: ", Rs," - System Transmittance: ",Ts, " - Sum: ", (Rs+Ts))
#print("p - polarized")
print("System Reflectance: ", Rp," - System Transmittance: ",Tp, " - Sum: ", (Rp+Tp))

