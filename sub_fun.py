import numpy as np

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
    sound_speed = calculate_seawater_sound_speed(temperature, salinity, depth, researcher)
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
    finf_values = calculate_elastic_sphere_form_function(
        frequency_factor, sound_speed, longitudinal_wave_speed, transverse_wave_speed, 
        medium_density, sphere_density, poisson_ratio)

    # Convert results to dB and display
    ts_values = 20 * np.log10(finf_values * diameter / 4000)
    #print('Ts = ', ts_values)
    return ts_values
    #print("vfa     sound_speed      finf")
    #print(frequency_factor, sound_speed, finf_values)

###########################################################

def rvf(value: str, default: str) -> str:
    """
    Return the default value if the input value is an empty string, otherwise return the input value.

    :param value: The input value to be checked. Expected to be a string.
    :param default: The default value to return if the input value is an empty string.
    :return: The input value if it is not an empty string, otherwise the default value.
    """
    return default if value == '' else value


def calculate_seawater_sound_speed(temperature: float = None, salinity: float = 34, depth: float = 100, researcher: str = 'mac') -> np.ndarray:
    """
    Calculate the speed of sound in seawater based on temperature, salinity, depth, 
    and the chosen empirical formula (Mackenzie or Medwin).

    :param temperature: Temperature in degrees Celsius. If None, it's requested via input.
    :param salinity: Salinity in practical salinity units (psu). Default is 34 psu.
    :param depth: Depth in meters. Default is 100 meters.
    :param researcher: Researcher's choice of formula: 'mac' for Mackenzie, 'med' for Medwin. Default is 'mac'.
    :return: Speed of sound in seawater calculated based on the given parameters. The return type is a numpy.ndarray.
    """
    # 如果温度temperature没有提供，通过输入获取temperature、salinity、depth和researcher的值
    if temperature is None:
        temperature = float(rvf(input('Temperature as [0:1:20] or 15(ret) in deg C ? '), 15))
        salinity = float(rvf(input('Salinity as [30:1:40] or 34(ret) in psu ? '), 34))
        depth = float(rvf(input('Depth as 50(ret) in m ? '), 50))
        researcher = rvf(input('Researcher: Mackenzie(mac,ret) or Medwin(med) ? '), 'mac')

    # 检查salinity是否为数组或列表，以支持向量化操作
    if isinstance(salinity, np.ndarray) or isinstance(salinity, list):
        # 确保temperature是列向量，以便于后续的向量化计算
        temperature = np.array(temperature).reshape(-1, 1)

        # 初始化列表sound_speeds，用于存储不同盐度下的声速
        sound_speeds = [temperature]
        for single_salinity in salinity:
            # 根据Mackenzie公式计算声速
            if researcher == 'mac':
                c = 1448.96 + 4.591*temperature - 0.05304*temperature**2 + \
                    2.374e-4*temperature**3 + 1.340*(single_salinity-35) + \
                    0.01630*depth + 1.675e-7*depth**2 - \
                    0.01025*temperature*(single_salinity-35) - 7.139e-13*temperature*depth**3
            # 根据Medwin公式计算声速
            elif researcher == 'med':
                c = 1449.2 + 4.6*temperature - 0.055*temperature**2 + 0.00029*temperature**3 + \
                    (1.34 - 0.010*temperature)*(single_salinity-35) + 0.016*depth
            # 将计算结果添加到列表sound_speeds中
            sound_speeds.append(c)
        # 将列表sound_speeds转换为numpy数组，并沿着列方向合并
        sound_speeds = np.concatenate(sound_speeds, axis=1)
    else:
        # 对于单个盐度值的情况，直接计算声速
        if researcher == 'mac':
            sound_speeds = 1448.96 + 4.591*temperature - 0.05304*temperature**2 + \
                2.374e-4*temperature**3 + 1.340*(salinity-35) + \
                0.01630*depth + 1.675e-7*depth**2 - \
                0.01025*temperature*(salinity-35) - 7.139e-13*temperature*depth**3
        elif researcher == 'med':
            sound_speeds = 1449.2 + 4.6*temperature - 0.055*temperature**2 + 0.00029*temperature**3 + \
                (1.34 - 0.010*temperature)*(salinity-35) + 0.016*depth

    # 返回计算出的声速
    return sound_speeds

################################################
# Error Handling Functions for Spherical Bessel Function Calculation
# These functions handle specific error scenarios encountered during the calculation
# of spherical Bessel functions. They output relevant error messages and return
# default values based on the error type.
# - `invalid` handles invalid arguments.
# - `not_accurate` handles scenarios where the calculation is not accurate.
# - `overflow` handles overflow errors.

def invalid(order, variable):
    #print("Argument of spherical_bessel_function is invalid. order={}, variable={}".format(order, variable))
    return 0

def not_accurate(order, variable):
    #print("Value of spherical_bessel_function is not accurate. order={}, variable={}".format(order, variable))
    return 0

def overflow(order, variable):
    #print("Value of spherical_bessel_function is overflow. order={}, variable={}".format(order, variable))
    return -1e30

################################################

def large1(order, variable):
    dx0 = 0.0000003
    if variable < 0.2:
        y = variable * variable
        w = 1 - y * (1 - 0.05 * y) / 6  # 近似 sin(variable) / variable
    else:
        m = np.floor(variable / np.pi)
        dx = variable - m * np.pi
        if dx < dx0:
            variable = m * np.pi + dx0  # 避免 integer * pi
        elif dx > np.pi - dx0:
            variable = (m + 1) * np.pi - dx0
        w = np.sin(variable) / variable

    if order < 0:
        return invalid(order, variable)
    elif order == 0:
        return w

    # 确定 l
    if variable >= 100:
        l = 0.02 * variable + 18
    elif 10 <= variable < 100:
        l = 0.1 * variable + 10
    elif 1 <= variable < 10:
        l = 0.5 * variable + 5
    else:
        l = 5

    nm = max(order, np.floor(variable)) + np.floor(l)
    z = 1 / variable
    t3 = 0
    t2 = 1e-35
    sj = 0  # 初始化 sj 以避免潜在的引用前未定义的错误
    for ii in range(1, int(nm) + 1):
        jj = nm - ii
        t1 = (jj + jj + 3) * z * t2 - t3
        if abs(order - jj) < 0.01:
            sj = t1
        if abs(t1) >= 1e25:
            t1 = t1 * 1e-25
            t2 = t2 * 1e-25
            sj = sj * 1e-25
        t3 = t2
        t2 = t1

    y = w / t1 * sj
    return y

################################################

def small1(order, variable):
    w = 1
    if order < 0:
        return invalid(order, variable)
    elif order == 0:
        return w
    elif 0 < order <= 10:
        t1 = 3
        t2 = 1
        for ii in range(1, order + 1):
            t3 = t2 * variable / t1
            t1 += 2
            t2 = t3
        return t3
    elif order > 10:
        return 0

################################################

def large2(order, variable):
    z = 1 / variable
    qn0 = -z * np.cos(variable)
    qn1 = z * (qn0 - np.sin(variable))

    if order < 0:
        return invalid(order, variable)
    elif order == 0:
        return qn0
    elif order == 1:
        return qn1
    elif order > 1:
        for ii in range(2, order + 1):
            qn2 = (ii + ii - 1) * z * qn1 - qn0
            if qn2 + 1e30 <= 0:
                return overflow(order, variable)
            qn0 = qn1
            qn1 = qn2
        return qn2
    
################################################
    
def small2(order, variable):
    z = 1 / variable
    if order < 0:
        return invalid(order, variable)
    elif order == 0:
        return -z
    elif order > 0:
        m = 30 / np.log10(z)
        if order - m + 1 > 0:
            return overflow(order, variable)
        qn0 = z
        qn1 = 1
        for ii in range(1, order + 1):
            qn2 = qn0 * qn1 * z
            qn1 += 2
            qn0 = qn2
        return -qn2
    
###########################################################

def spherical_bessel_function(kind, order, variable):
    """
    Calculate spherical Bessel functions and their first and second derivatives.
    
    :param kind: Type of the function, 1 for first kind (j), 2 for second kind (n).
    :param order: Order of the function.
    :param variable: The variable to evaluate the function at.
    :return: Tuple of (function value, first derivative, second derivative).
    """
    # Check for invalid inputs
    if order < 0 or variable < 0:
        raise ValueError("Invalid input: order and variable must be non-negative")

    # Initialize arrays to store values for the function and its derivatives
    max_order = order + 2  # Calculate up to two orders higher for derivatives
    function_value = np.zeros(max_order + 1)
    first_derivative = np.zeros(max_order + 1)
    second_derivative = np.zeros(max_order + 1)

    for current_order in range(max_order + 1):
        if kind == 1:  # first kind
            if current_order > 30000:
                return not_accurate(current_order, variable)
            if variable < 0:
                return invalid(current_order, variable)
            elif 0 <= variable < 7e-4:
                # function_value = small1(order, variable)
                function_value[current_order] = small1(current_order, variable)
            elif 7e-4 <= variable < 3e4:
                # function_value = large1(order, variable)
                function_value[current_order] = large1(current_order, variable)
            else:
                return not_accurate(current_order, variable)

        elif kind == 2:  # second kind
            if current_order >= 3e4:
                return not_accurate(current_order, variable)
            if variable < 0:
                return invalid(current_order, variable)
            elif variable == 0:
                function_value[current_order] = -1e30
            elif 0 < variable <= 7e-4:
                # function_value = small2(order, variable)
                function_value[current_order] = small2(current_order, variable)
            else:
                # function_value = large2(order, variable)
                function_value[current_order] = large2(current_order, variable)
    '''

    # Calculate derivatives if required
    function_value_plus_one = spherical_bessel_function(kind, order + 1, variable)
    first_derivative = order / variable * function_value - function_value_plus_one

    function_value_plus_two = spherical_bessel_function(kind, order + 2, variable)
    second_derivative = ((2 * order + 1) / variable * first_derivative -
                            order * (order + 2) / (variable ** 2) * function_value +
                            function_value_plus_two)

    return function_value, first_derivative, second_derivative
    '''


    # Calculate derivatives
    for i in range(max_order):
        first_derivative[i] = i / variable * function_value[i] - function_value[i + 1]

    for i in range(max_order - 1):
        second_derivative[i] = ((2 * i + 1) / variable * first_derivative[i] -
                                 i * (i + 2) / (variable ** 2) * function_value[i] +
                                 function_value[i + 2])

    # Return the values for the requested order
    return function_value[order], first_derivative[order], second_derivative[order]

###########################################################

def calculate_elastic_sphere_form_function(frequency=580, medium_speed=1500, longitudinal_wave_speed=6867, transverse_wave_speed=4161.2, medium_density=1027, sphere_density=14900, poisson_ratio=0.2095):
    """
    Calculate the form function of an elastic sphere.

    :param frequency: Frequency in Hz.
    :param medium_speed: Speed of sound in the surrounding medium (like water) in m/s.
    :param longitudinal_wave_speed: Longitudinal wave speed in the sphere material in m/s.
    :param transverse_wave_speed: Transverse wave speed in the sphere material in m/s.
    :param medium_density: Density of the surrounding medium in kg/m^3.
    :param sphere_density: Density of the sphere material in kg/m^3.
    :param poisson_ratio: Poisson's ratio of the sphere material.
    :return: The absolute value of the form function of the elastic sphere.

    The calculation is based on the theory by Faran as referenced in Miyanohana et al. (1993).
    """

    wave_number = 2 * np.pi * frequency / medium_speed
    nmax = int(np.round(3 * wave_number))
    longitudinal_wave_number = 2 * np.pi * frequency / longitudinal_wave_speed
    transverse_wave_number = 2 * np.pi * frequency / transverse_wave_speed

    amplitude_coefficients = np.zeros(nmax + 1, dtype=np.complex128)

    for mode in range(nmax + 1):
        jx, jxd, _ = spherical_bessel_function(1, mode, wave_number)
        nx, nxd, _ = spherical_bessel_function(2, mode, wave_number)
        jx1, jx1d, jx1dd = spherical_bessel_function(1, mode, longitudinal_wave_number)
        jx2, jx2d, jx2dd = spherical_bessel_function(1, mode, transverse_wave_number)

        nn = mode * (mode + 1)
        f1 = longitudinal_wave_number * jx1d - jx1
        f2 = (nn - 2) * jx2 + transverse_wave_number ** 2 * jx2dd
        f3 = longitudinal_wave_number * jx1d / f1
        f4 = 2 * nn * jx2 / f2
        f5 = longitudinal_wave_number ** 2 * (poisson_ratio / (1 - 2 * poisson_ratio) * jx1 - jx1dd) / f1
        f6 = 2 * nn * (jx2 - transverse_wave_number * jx2d) / f2
        fn = medium_density / sphere_density * transverse_wave_number ** 2 / 2 * (f3 - f4) / (f5 - f6)

        etan = np.arctan(-(jx * fn - wave_number * jxd) / (nx * fn - wave_number * nxd))
        amplitude_coefficients[mode] = (-1) ** mode * (2 * mode + 1) * np.sin(etan) * np.exp(-1j * etan)

    form_function = 2 / wave_number * np.abs(np.sum(amplitude_coefficients))
    return form_function






