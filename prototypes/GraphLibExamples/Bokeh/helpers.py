import numpy as np
from typing import List,Tuple

def veclist_to_coordlists(veclist):
    xs=[v[0] for v in veclist]
    ys=[v[1] for v in veclist]
    return xs,ys
    #return list(zip(*veclist))

def norm(vec):
    return np.sqrt(np.dot(vec,vec))

def with_arrow_head(
        xs: List[float],
        ys: List[float], 
        width: float = 0.5,
        length: float =0.5,
    )->Tuple[List[float],List[float]]:
    
    last_point=np.array([
        xs[-1],
        ys[-1]
    ])
    second_last_point=np.array([
        xs[-2],
        ys[-2]
    ])
    last_vector=last_point-second_last_point
    last_vector_n=last_vector/norm(last_vector)
    print('last_vector',last_vector)
    left_head_vector=np.array([
        -last_vector_n[1],
         last_vector_n[0]
    ])    
    print('left_head_vector',left_head_vector)
    right_head_vector = -left_head_vector
    left_head_point=last_point -length*last_vector_n +width*left_head_vector
    
    right_head_point = last_point - length*last_vector_n  +  width*right_head_vector
    #we go first from the tip to left head
    left_stroke=[left_head_point,last_point]
    right_stroke=[right_head_point,last_point]

    
    left_xs,left_ys = veclist_to_coordlists(left_stroke)
    right_xs,right_ys = veclist_to_coordlists(right_stroke)
    xss=xs+left_xs+right_xs
    yss=ys+left_ys+right_ys
    return xss,yss
