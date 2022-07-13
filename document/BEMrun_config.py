stl_unit = 1e-3 # meters

# this is the unit for our Mesh structure
# in the following code, this is the default unit for python variables
mesh_unit = 1e-3 # meters
filename = 'glasstrapcolor.stl'


SAVE_PLOTS = False
SHOW_PLOTS = True
USE_MULTIPROCESSING = False

# assign a name for each color
electrodes = {
    ('bem1', 'DC1'),
    ('bem2', 'DC2'),
    ('bem3', 'DC3'),
    ('bem4', 'DC4'),
    ('bem5', 'DC5'),
    ('bem6', 'DC6'),
    ('bem7', 'DC7'),
    ('bem8', 'DC8'),
    ('bem9', 'DC9'),
    ('bem10', 'DC10'),
    ('bem11', 'DC11'),
    ('bem12', 'DC12'),
    ('bem13', 'DC13'),
    ('bem14', 'DC14'),
    ('bem15', 'DC15'),
    ('bem16', 'DC16'),
    ('bem17', 'DC17'),
    ('bem18', 'DC18'),
    ('bem19', 'DC19'),
    ('bem20', 'DC20'),
    ('bem21', 'DC21'),
    ('bem22', 'RF'),
    ('bem23', 'GND'),
}

xl = 0*1e-3
yl = 10*1e-3
zl = 2 + 1*75*1e-3
