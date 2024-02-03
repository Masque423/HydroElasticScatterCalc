import tkinter as tk
from tkinter import ttk
import sub_fun

def calculate_ts(event=None):
    try:
        # 获取下拉列表及输入框中的各参数值
        material = str(material_combobox.get())
        temperature = float(temperature_entry.get())
        salinity = float(salinity_entry.get())
        depth = float(depth_entry.get())
        frequency = float(frequency_entry.get())
        diameter = float(diameter_entry.get())
        researcher = str(researcher_combobox.get())
        
        # 计算散射特性
        TS = sub_fun.hyfro_elastic_scatter_calc(material, temperature, salinity, depth, frequency, diameter, researcher)
        # 在结果标签中显示总和
        result_label.config(text="Ts: " + str(TS))
    except ValueError:
        # 如果输入不是有效的整数，显示错误信息
        result_label.config(text="Please enter valid integers.")


# 定义一个函数用于居中窗口
def center_window(width=300, height=250):
    # 获取屏幕宽度和高度
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # 计算窗口的 x 和 y 坐标
    x = (screen_width / 2) - (width / 2)
    y = (screen_height / 2) - (height / 2)

    # 设置窗口的几何尺寸（宽、高、x 坐标、y 坐标）
    root.geometry('%dx%d+%d+%d' % (width, height, x, y))
    
# 创建主窗口实例
root = tk.Tk()
root.title("ScatterCalc")

# 将窗口置于屏幕中央
center_window(300, 500)

# 创建 material 输入的标签，并显示
material_label = ttk.Label(root, text="Material (Cu, WC, Steel, Nylon, etc.): ")
material_label.pack()
# 创建 material 的下拉列表
material_options = ['Cu', 'WC', 'Steel', 'Nylon']
material_combobox = ttk.Combobox(root, values=material_options)
material_combobox.pack()
material_combobox.bind("<Return>", calculate_ts)
material_combobox.set(material_options[0])  # 设置默认值为列表中的第一个值


# 创建 temperature 输入的标签，并显示
temperature_label = ttk.Label(root, text = "Water temperature [degree C]: ")
temperature_label.pack()
# 创建 temperature 的输入框并设置为窗口打开时的焦点
temperature_entry = tk.Entry(root)
temperature_entry.pack()
temperature_entry.focus_set()
# 绑定回车键事件，当按下回车时调用 calculate_ts 函数
temperature_entry.bind("<Return>", calculate_ts)


# 创建 salinity 输入的标签和输入框，并显示
salinity_label = ttk.Label(root, text = "Salinity [psu]: ")
salinity_label.pack()
salinity_entry = tk.Entry(root)
salinity_entry.pack()
salinity_entry.bind("<Return>", calculate_ts)

# 创建 depth 输入的标签和输入框，并显示
depth_label = ttk.Label(root, text = "Depth [m]: ")
depth_label.pack()
depth_entry = tk.Entry(root)
depth_entry.pack()
depth_entry.bind("<Return>", calculate_ts)


# 创建 frequency 输入的标签和输入框，并显示
frequency_label = ttk.Label(root, text = "Frequency [kHz]: ")
frequency_label.pack()
frequency_entry = tk.Entry(root)
frequency_entry.pack()
frequency_entry.bind("<Return>", calculate_ts)

# 创建 diameter 输入的标签和输入框，并显示
diameter_label = ttk.Label(root, text = "Diameter [mm]: ")
diameter_label.pack()
diameter_entry = tk.Entry(root)
diameter_entry.pack()
diameter_entry.bind("<Return>", calculate_ts)


# 创建 researcher 输入的标签，并显示
researcher_label = ttk.Label(root, text = "researcher: ")
researcher_label.pack()
# 创建 researcher 的下拉列表
researcher_option = ['mac', 'WC']
researcher_combobox = ttk.Combobox(root, values=researcher_option)
researcher_combobox.pack()
researcher_combobox.bind("<Return>", calculate_ts)
researcher_combobox.set(researcher_option[0])  # 设置默认值为列表中的第一个值


# 创建一个标签，用于显示计算结果
result_label = tk.Label(root, text="TS: ")
result_label.pack()


# result_entry = tk.Entry(root)
# esult_entry.pack()
# result_entry.bind("<Return>", calculate_ts)


# 启动 Tkinter 事件循环
root.mainloop()
