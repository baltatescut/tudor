import pathlib


def load_system(file_path: pathlib.Path) -> tuple[list[list[float]], list[float]]:
    matrix_a: list[list[float]] = []
    vector_b: list[float] = []

    with open(file_path, "r") as file_obj:
        for line in file_obj:
            if "=" not in line:
                continue
            left_part, right_part = line.split("=")
            vector_b.append(float(right_part))

            current_row = []
            sign_val = 1.0
            for token in left_part.split():
                if token == "+":
                    sign_val = 1.0
                elif token == "-":
                    sign_val = -1.0
                else:
                    val_str = token[:-1]
                    if val_str == "":
                        parsed_val = 1.0
                    else:
                        parsed_val = float(val_str)
                    current_row.append(parsed_val * sign_val)
            matrix_a.append(current_row)
    return matrix_a, vector_b


def determinant(mat: list[list[float]]) -> float:
    a11, a12, a13 = mat[0]
    a21, a22, a23 = mat[1]
    a31, a32, a33 = mat[2]
    det_val = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31)
    return float(det_val)


def trace(mat: list[list[float]]) -> float:
    trace_val = mat[0][0] + mat[1][1] + mat[2][2]
    return float(trace_val)


def norm(vec: list[float]) -> float:
    sum_squares = vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2
    return float(sum_squares ** 0.5)


def transpose(mat: list[list[float]]) -> list[list[float]]:
    size = len(mat)
    result_mat = []
    for i in range(size):
        new_row = []
        for j in range(size):
            new_row.append(mat[j][i])
        result_mat.append(new_row)
    return result_mat


def multiply(mat: list[list[float]], vec: list[float]) -> list[float]:
    res_vec = []
    size = len(mat)
    for i in range(size):
        val = mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2]
        res_vec.append(val)
    return res_vec


def solve_cramer(mat: list[list[float]], vec: list[float]) -> list[float]:
    det_a = determinant(mat)

    mat_x = [[vec[0], mat[0][1], mat[0][2]],
             [vec[1], mat[1][1], mat[1][2]],
             [vec[2], mat[2][1], mat[2][2]]]

    mat_y = [[mat[0][0], vec[0], mat[0][2]],
             [mat[1][0], vec[1], mat[1][2]],
             [mat[2][0], vec[2], mat[2][2]]]

    mat_z = [[mat[0][0], mat[0][1], vec[0]],
             [mat[1][0], mat[1][1], vec[1]],
             [mat[2][0], mat[2][1], vec[2]]]

    x = determinant(mat_x) / det_a
    y = determinant(mat_y) / det_a
    z = determinant(mat_z) / det_a
    return [x, y, z]


def minor(mat: list[list[float]], row_idx: int, col_idx: int) -> list[list[float]]:
    min_matrix = []
    size = len(mat)
    for r in range(size):
        if r == row_idx:
            continue
        row_items = []
        for c in range(size):
            if c == col_idx:
                continue
            row_items.append(mat[r][c])
        min_matrix.append(row_items)
    return min_matrix


def cofactor(mat: list[list[float]]) -> list[list[float]]:
    size = len(mat)
    cof_matrix = []
    for i in range(size):
        row_items = []
        for j in range(size):
            sub_mat = minor(mat, i, j)
            det_sub = sub_mat[0][0] * sub_mat[1][1] - sub_mat[0][1] * sub_mat[1][0]
            row_items.append(((-1) ** (i + j)) * det_sub)
        cof_matrix.append(row_items)
    return cof_matrix


def adjoint(mat: list[list[float]]) -> list[list[float]]:
    cof_mat = cofactor(mat)
    adj_mat = transpose(cof_mat)
    return adj_mat


def solve(mat: list[list[float]], vec: list[float]) -> list[float]:
    det_a = determinant(mat)
    adj_a = adjoint(mat)
    size = len(mat)
    inv_mat = []

    for i in range(size):
        row_items = []
        for j in range(size):
            row_items.append((1 / det_a) * adj_a[i][j])
        inv_mat.append(row_items)

    solution_vec = multiply(inv_mat, vec)
    return solution_vec