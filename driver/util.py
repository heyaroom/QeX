import re

def name_to_num(name):
    num = int(re.sub("\\D", "", name))
    return num
def num_to_alpha(num):
    alpha = chr(num+64)
    return alpha
def name_to_alpha(name):
    num = name_to_num(name)
    alpha = num_to_alpha(num)
    return alpha