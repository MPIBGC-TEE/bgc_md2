from pathlib import Path
from functools import reduce

def rp(s,t):
    """produces the relative path from source to target
    as in an arguemt for chdir"""
    sp=s.absolute().parts
    tp=t.absolute().parts
    common, rest_s, rest_t = split_3(sp,tp)
    return Path().joinpath(
        *(
            [".." for el in rest_s]
            +[el for el in rest_t]
        )
    )

def split_3(sp,tp):
    """looking at two paths returns 
    1. the common part, 
    2. the src path without the common part and 
    3. the target part without the common part
    """
    return split_3_rec((),sp,tp)

def split_3_rec(common,sp,tp):
    if min(len(sp),len(tp)) == 0:
        return common,sp,tp
    elif sp[0]!=tp[0]:
        return common,sp,tp
    else:
        return split_3_rec(
            common+(sp[0],),
            sp[1:],
            tp[1:]
        )

