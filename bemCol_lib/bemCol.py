#convert rgb value in fusion360 to color index in stl
def _fusion360_RGBA32_to_stlRGB16(r:int,g:int,b:int):
    if r<0 or r>255 or g<0 or g>255 or b<0 or b>255:
        raise ValueError('not a rgb value')
    def _8bit_to_5bit(x):
        if x == 255:
            return 31
        else:
            return int(x/255*31)
    return _8bit_to_5bit(b)*32**2 + _8bit_to_5bit(g)*32 + _8bit_to_5bit(r)

rgb2stl = _fusion360_RGBA32_to_stlRGB16

# standart bemCol rgb values
_bemCol_rgb = {'bem1': (160, 16, 112), 'bem2': (128, 160, 16), 'bem3': (64, 128, 224), 'bem4': (128, 224, 80), 'bem5': (176, 208, 48), 'bem6': (16, 144, 128), 'bem7': (96, 32, 64), 'bem8': (192, 96, 16), 'bem9': (112, 224, 48), 'bem10': (160, 96, 112), 'bem11': (112, 80, 224), 'bem12': (176, 80, 128), 'bem13': (48, 176, 80), 'bem14': (176, 96, 48), 'bem15': (128, 208, 224), 'bem16': (64, 128, 160), 'bem17': (112, 16, 160), 'bem18': (16, 64, 48), 'bem19': (64, 240, 80), 'bem20': (144, 96, 208), 'bem21': (240, 64, 96), 'bem22': (176, 144, 128), 'bem23': (112, 128, 32), 'bem24': (176, 160, 96), 'bem25': (224, 144, 48), 'bem26': (224, 64, 112), 'bem27': (64, 16, 16), 'bem28': (128, 80, 224), 'bem29': (160, 176, 32), 'bem30': (144, 240, 224), 'bem31': (32, 48, 112), 'bem32': (160, 96, 240), 'bem33': (160, 48, 32), 'bem34': (208, 208, 32), 'bem35': (128, 48, 208), 'bem36': (64, 96, 112), 'bem37': (112, 128, 64), 'bem38': (176, 80, 48), 'bem39': (160, 144, 48), 'bem40': (224, 32, 80)}

# a dict {'bemxx':stlColorValue}
bemCol = dict([(_key,rgb2stl(_bemCol_rgb[_key][0],_bemCol_rgb[_key][1],_bemCol_rgb[_key][2])) for _key in _bemCol_rgb.keys()])

'''
or 
bemCol = {'bem1': 13363, 'bem2': 1647, 'bem3': 28135, 'bem4': 10095, 'bem5': 5941, 'bem6': 15905, 'bem7': 7275, 'bem8': 1399, 'bem9': 5997, 'bem10': 13683, 'bem11': 27949, 'bem12': 15669, 'bem13': 9893, 'bem14': 5493, 'bem15': 28463, 'bem16': 19943, 'bem17': 19501, 'bem18': 5345, 'bem19': 10151, 'bem20': 25969, 'bem21': 11517, 'bem22': 15925, 'bem23': 3565, 'bem24': 11893, 'bem25': 5691, 'bem26': 13563, 'bem27': 1063, 'bem28': 27951, 'bem29': 3763, 'bem30': 28593, 'bem31': 13475, 'bem32': 30067, 'bem33': 3251, 'bem34': 3897, 'bem35': 25775, 'bem36': 13671, 'bem37': 7661, 'bem38': 5429, 'bem39': 5683, 'bem40': 9339}
'''

_t = 2**5
meshlab_bemCol = dict(zip(bemCol.keys(),[_t**3 + (_t**2)*(k%_t) + _t*((k//_t)%_t) + k//(_t**2) for k in bemCol.values()]))

print('meshlab_bemCol',meshlab_bemCol)


