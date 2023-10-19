import cv2
import numpy as np
import scipy.io as scio
# 4邻域的连通域和 8邻域的连通域
# [row, col]
NEIGHBOR_HOODS_4 = True
OFFSETS_4 = [[0, -1], [-1, 0], [0, 0], [1, 0], [0, 1]]

NEIGHBOR_HOODS_8 = False
OFFSETS_8 = [[-1, -1], [0, -1], [1, -1],
             [-1,  0], [0,  0], [1,  0],
             [-1,  1], [0,  1], [1,  1]]



def reorganize(binary_img: np.array):
    index_map = []
    points = []
    index = -1
    rows, cols = binary_img.shape
    for row in range(rows):
        for col in range(cols):
            var = binary_img[row][col]
            if var < 0.5:
                continue
            if var in index_map:
                index = index_map.index(var)
                num = index + 1
            else:
                index = len(index_map)
                num = index + 1
                index_map.append(var)
                points.append([])
            binary_img[row][col] = num
            points[index].append([row, col])
    return binary_img, points



def neighbor_value(binary_img: np.array, offsets, reverse=False):
    rows, cols = binary_img.shape
    label_idx = 0
    rows_ = [0, rows, 1] if reverse == False else [rows-1, -1, -1]
    cols_ = [0, cols, 1] if reverse == False else [cols-1, -1, -1]
    for row in range(rows_[0], rows_[1], rows_[2]):
        for col in range(cols_[0], cols_[1], cols_[2]):
            label = 256
            if binary_img[row][col] < 0.5:
                continue
            for offset in offsets:
                neighbor_row = min(max(0, row+offset[0]), rows-1)
                neighbor_col = min(max(0, col+offset[1]), cols-1)
                neighbor_val = binary_img[neighbor_row, neighbor_col]
                if neighbor_val < 0.5:
                    continue
                label = neighbor_val if neighbor_val < label else label
            if label == 255:
                label_idx += 1
                label = label_idx
            binary_img[row][col] = label
    return binary_img

# binary_img: bg-0, object-255; int
def Two_Pass(binary_img: np.array, neighbor_hoods):
    if neighbor_hoods == 'NEIGHBOR_HOODS_4':
        offsets = [[0, -1], [-1, 0], [0, 0], [1, 0], [0, 1]]
    elif neighbor_hoods == 'NEIGHBOR_HOODS_8':
        offsets = [[-1, -1], [0, -1], [1, -1],
             [-1,  0], [0,  0], [1,  0],
             [-1,  1], [0,  1], [1,  1]]
    else:
        raise ValueError

    binary_img = neighbor_value(binary_img, offsets, False)
    binary_img = neighbor_value(binary_img, offsets, True)

    return binary_img



def recursive_seed(binary_img: np.array, seed_row, seed_col, offsets, num, max_num=100):
    rows, cols = binary_img.shape
    binary_img[seed_row][seed_col] = num
    for offset in offsets:
        neighbor_row = min(max(0, seed_row+offset[0]), rows-1)
        neighbor_col = min(max(0, seed_col+offset[1]), cols-1)
        var = binary_img[neighbor_row][neighbor_col]
        if var < max_num:
            continue
        binary_img = recursive_seed(binary_img, neighbor_row, neighbor_col, offsets, num, max_num)
    return binary_img

# max_num 表示连通域最多存在的个数
def Seed_Filling(binary_img, neighbor_hoods, max_num=1000000):
    if neighbor_hoods == 'NEIGHBOR_HOODS_4':
        offsets = [[0, -1], [-1, 0], [0, 0], [1, 0], [0, 1]]
    elif neighbor_hoods == 'NEIGHBOR_HOODS_8':
        offsets = [[-1, -1], [0, -1], [1, -1],
             [-1,  0], [0,  0], [1,  0],
             [-1,  1], [0,  1], [1,  1]]
    else:
        raise ValueError

    num = 1
    rows, cols = binary_img.shape
    for row in range(rows):
        for col in range(cols):
            var = binary_img[row][col]
            if var <= max_num:
                continue
            binary_img = recursive_seed(binary_img, row, col, offsets, num, max_num=1000000)
            num += 1
    return binary_img


def cal_pccs(x, y, n):
    """
    warning: data format must be narray
    :param x: Variable 1
    :param y: The variable 2
    :param n: The number of elements in x
    :return: pccs
    """
    sum_xy = np.sum(np.sum(x*y))
    sum_x = np.sum(np.sum(x))
    sum_y = np.sum(np.sum(y))
    sum_x2 = np.sum(np.sum(x*x))
    sum_y2 = np.sum(np.sum(y*y))
    pcc = (n*sum_xy-sum_x*sum_y)/np.sqrt((n*sum_x2-sum_x*sum_x)*(n*sum_y2-sum_y*sum_y))
    return pcc


def four_cc_label(img):
    height, width = img.shape
    label = np.zeros((height, width), dtype=np.int32)
    LUT = np.zeros(height * width,dtype=np.uint8)

    COLORS = range(0,1000000)
    out = np.zeros((height, width), dtype=np.uint8)
    label[img[:,:] > 0] = 1

    n = 1
    for y in range(height):
        for x in range(width):
            if y%500==0:
                if x%500==0:
                    print(y,' --- ',height,' --- ',x,' --- ',width)
            if label[y, x] == 0:
                continue
            c2 = label[max(y - 1, 0), min(x + 1, width - 1)]
            c3 = label[max(y - 1, 0), x]
            c4 = label[max(y - 1, 0), max(x - 1, 0)]
            c5 = label[y, max(x - 1, 0)]
            if c3 < 2 and c5 < 2 and c2 < 2 and c4 < 2:
                n += 1
                label[y, x] = n
            else:
                _vs = [c3, c5, c2, c4]
                vs = [a for a in _vs if a > 1]
                v = min(vs)
                label[y, x] = v

                minv = v
                for _v in vs:
                    if LUT[_v] != 0:
                        minv = min(minv, LUT[_v])
                for _v in vs:
                    LUT[_v] = minv

    count = 1
    for l in range(2, n + 1):
        flag = True
        for i in range(n + 1):
            if LUT[i] == l:
                if flag:
                    count += 1
                    flag = False
                LUT[i] = count


    for i, lut in enumerate(LUT[2:]):
        if i%500==0:
            print('LUT',i,len(LUT[2:]))
        out[label == (i + 2)] = COLORS[lut - 2]
    return out



if __name__ == "__main__":
    binary_img = np.zeros((4, 7), dtype=np.int16)
    index = [[0, 2], [0, 5],
            [1, 0], [1, 1], [1, 2], [1, 4], [1, 5], [1, 6],
            [2, 2], [2, 5],
            [3, 1], [3, 2], [3, 4], [3, 6]]
    for i in index:
        binary_img[i[0], i[1]] = np.int16(255)

    print("原始二值图像")
    print(binary_img)

    print("Two_Pass")
    binary_img = Two_Pass(binary_img, NEIGHBOR_HOODS_8)
    binary_img, points = reorganize(binary_img)
    print(binary_img, points)

    print("Seed_Filling")
    binary_img = Seed_Filling(binary_img, NEIGHBOR_HOODS_8)
    binary_img, points = reorganize(binary_img)
    print(binary_img, points)