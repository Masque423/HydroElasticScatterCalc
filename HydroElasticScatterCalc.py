"""
scatesf: Scattering characteristics of elatic sphere
计算并显示弹性球体在水中的散射特性
"""
import numpy as np
import sub_fun

def hyfro_elastic_scatter_calc(material:str, temperature:float, salinity:float, depth:float, frequency:float, diameter:float, researcher:str):

    '''
    material = input('Material (Cu, WC, Steel, Nylon, etc.): ')
    temperature = float(input('Water temperature [degree C]: '))
    salinity = float(input('Salinity [psu]: '))
    depth = float(input('Depth [m]: '))
    frequency = float(input('Frequency [kHz]: '))
    diameter = float(input('Diameter [mm]: '))
    researcher = 'mac'

    material = 'Cu'
    temperature = 10
    salinity = 34
    depth = 100
    frequency = 20
    diameter = 10
    researcher = 'mac'
    '''
    # 计算水中的声速
    sound_speed = sub_fun.calculate_seawater_sound_speed(temperature, salinity, depth, researcher)
    #c_values = [sound_speed, 1, sound_speed]

    frequency_factor = frequency * diameter / 2
    #frequency_array = [frequency_factor, 1, frequency_factor]

    if material == 'WC':
        longitudinal_wave_speed = 6867
        transverse_wave_speed = 4161.2
        medium_density = 1000
        sphere_density = 14900
        poisson_ratio = 0.2095
    elif material == 'Cu':
        longitudinal_wave_speed = 4760
        transverse_wave_speed = 2288.5
        medium_density = 1000
        sphere_density = 8947
        poisson_ratio = 0.3497
    elif material == 'Steel':
        longitudinal_wave_speed = 5764.5
        transverse_wave_speed = 3230.8
        medium_density = 1000
        sphere_density = 7830
        poisson_ratio = 0.271
    elif material == 'Nylon':
        longitudinal_wave_speed = 2620
        transverse_wave_speed = 1070
        medium_density = 1024
        sphere_density = 1110
        poisson_ratio = 0.4
    else:
        raise ValueError("Unknown material")

    # Calculate scattering effect
    #finf_values = np.zeros((len(frequency_array), len(c_values)))
    #finf_values = np.zeros(1, 1)
    #for jj, vc in enumerate(c_values):
    #    jj += 1
    #    for ii, vfa in enumerate(frequency_array):
    #        ii += 1
    #vc = sound_speed
    #vfa = frequency_factor
    finf_values = sub_fun.calculate_elastic_sphere_form_function(
        frequency_factor, sound_speed, longitudinal_wave_speed, transverse_wave_speed, 
        medium_density, sphere_density, poisson_ratio)

    # Convert results to dB and display
    ts_values = 20 * np.log10(finf_values * diameter / 4000)
    print('Ts = ', ts_values)
    print("vfa     sound_speed      finf")
    print(frequency_factor, sound_speed, finf_values)


