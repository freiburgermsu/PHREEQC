import subprocess
import rosspy
import pandas
import os, re, time

simulation_folder = r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\PHREEQC\PHREEQC Test\2022-01-4-ROSSpy-red_sea-transport-pitzer-scaling-all_distance-LinPerm"
if not os.path.exists(simulation_folder):
    os.mkdir(r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\PHREEQC\PHREEQC Test\2022-01-4-ROSSpy-red_sea-transport-pitzer-scaling-all_distance-LinPerm")

bat_path = r"phreeqc.bat"

# input_path = r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\PHREEQC\PHREEQC Test\input.pqi"
input_path = input('what is the input path?')

# output_path = r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\PHREEQC\PHREEQC Test\2022-01-04-ROSSpy--transport-pitzer-scaling-all_distance-LinPerm\output.pqo"
output_path = input('what is the output path?')
max_argument_length = len(r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\ROSS\examples\scaling\scale_validation\2021-12-18-ROSSpy-red_sea-transport-pitzer-scaling-all_distance-LinPerm\2021-12-28-ROSSpy--transport-pitzer-scaling-all_")
if len(output_path) > max_argument_length:
    excessive_characters = output_path[max_argument_length:]
    output_path = re.sub('(input.pqi)', 'output.pqo', input_path)
    print(f'''\n--> ERROR: The output file path was abridged to {output_path} to maintain validity as an argument for the batch PHREEQC software.\n''')

# database_path = r"C:\Users\Andrew Freiburger\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.8_qbz5n2kfra8p0\LocalCache\local-packages\Python38\site-packages\rosspy\databases\pitzer.dat"
database_path = input('what is the database path?')

proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
command = str.encode(bat_path + " \"" + input_path + "\" \"" + output_path + "\" \"" + database_path + "\"\n") 
proc.stdin.write(command)
#time.sleep(4)
#proc.stdin.write(str.encode(database_path))
proc.stdin.close()  
proc.wait()