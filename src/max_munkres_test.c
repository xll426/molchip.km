
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>    // 使用 fabs 函数
#include <stdbool.h>

// 定义最大矩阵大小
#define MAX_SIZE 100

// 定义一个特殊的值来表示DISALLOWED
#define DISALLOWED_VAL DBL_MAX

// 定义标记类型
#define UNMARKED 0
#define STARRED 1
#define PRIMED 2

// 定义一个结构体来存储Munkres算法的状态
typedef struct {
    float C[MAX_SIZE][MAX_SIZE];           // 成本矩阵（会被修改）
    float original_C[MAX_SIZE][MAX_SIZE];  // 原始成本矩阵（保留不变）
    int marked[MAX_SIZE][MAX_SIZE];        // 标记矩阵（STARRED, PRIMED）
    bool row_covered[MAX_SIZE];            // 行覆盖标记
    bool col_covered[MAX_SIZE];            // 列覆盖标记
    int path[MAX_SIZE * 2][2];             // 路径矩阵
    int n;                                  // 矩阵大小
    int Z0_r;                               // 路径起始点行
    int Z0_c;                               // 路径起始点列
    float lx[MAX_SIZE];                    // 行标签
    float ly[MAX_SIZE];                    // 列标签
} Munkres;

// 定义一个结构体来存储结果
typedef struct {
    int row;
    int col;
} Assignment;

// 定义一个结构体来存储测试用例
typedef struct {
    float matrix[MAX_SIZE][MAX_SIZE];
    int rows;
    int cols;
    float expected_cost;
} TestCase;

// 打印矩阵的函数
void print_matrix(float matrix[MAX_SIZE][MAX_SIZE], int rows, int cols, const char* msg) {
    if (msg != NULL) {
        printf("%s\n", msg);
    }
    for (int i = 0; i < rows; i++) {
        printf("[");
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] == DISALLOWED_VAL) {
                printf("D"); // D 表示 DISALLOWED
            } else {
                printf("%.4lf", matrix[i][j]);
            }
            if (j < cols - 1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
    printf("\n");
}

// 数值取反功能
void negate_matrix(float matrix[MAX_SIZE][MAX_SIZE], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] != DISALLOWED_VAL) {
                matrix[i][j] = -matrix[i][j];
            }
        }
    }
}

// 创建一个新的矩阵并填充它，同时保留原始矩阵
void pad_matrix(Munkres* munkres, float input_matrix[][MAX_SIZE], int input_rows, int input_cols) {
    // 找到最大行和列
    int max_dim = input_rows > input_cols ? input_rows : input_cols;
    munkres->n = max_dim;

    // 填充矩阵，使用0.0作为填充值，并保留原始矩阵
    for (int i = 0; i < munkres->n; i++) {
        for (int j = 0; j < munkres->n; j++) {
            if (i < input_rows && j < input_cols) {
                munkres->C[i][j] = input_matrix[i][j];
                munkres->original_C[i][j] = input_matrix[i][j];
            } else {
                munkres->C[i][j] = 0.0; // 填充值为0.0
                munkres->original_C[i][j] = DISALLOWED_VAL; // 保留原始矩阵为DISALLOWED_VAL
            }
        }
    }
}

// 初始化Munkres结构体
void initialize(Munkres* munkres) {
    for (int i = 0; i < munkres->n; i++) {
        munkres->row_covered[i] = false;
        munkres->col_covered[i] = false;
        munkres->lx[i] = 0.0;
        munkres->ly[i] = 0.0;
        for (int j = 0; j < munkres->n; j++) {
            munkres->marked[i][j] = UNMARKED;
        }
    }
    memset(munkres->path, 0, sizeof(munkres->path));
    munkres->Z0_r = 0;
    munkres->Z0_c = 0;
}

// 查找未覆盖的零
bool find_a_zero(Munkres* munkres, int* row, int* col) {
    float epsilon = 1e-6;
    for (int i = 0; i < munkres->n; i++) {
        for (int j = 0; j < munkres->n; j++) {
            if (fabs(munkres->C[i][j]) < epsilon && !munkres->row_covered[i] && !munkres->col_covered[j] && munkres->original_C[i][j] != DISALLOWED_VAL) {
                *row = i;
                *col = j;
                return true;
            }
        }
    }
    return false;
}

// 查找行中第一个星号零
int find_star_in_row(Munkres* munkres, int row) {
    for (int j = 0; j < munkres->n; j++) {
        if (munkres->marked[row][j] == STARRED) {
            return j;
        }
    }
    return -1;
}

// 查找列中第一个星号零
int find_star_in_col(Munkres* munkres, int col) {
    for (int i = 0; i < munkres->n; i++) {
        if (munkres->marked[i][col] == STARRED) {
            return i;
        }
    }
    return -1;
}

// 查找行中第一个标记零
int find_prime_in_row(Munkres* munkres, int row) {
    for (int j = 0; j < munkres->n; j++) {
        if (munkres->marked[row][j] == PRIMED) {
            return j;
        }
    }
    return -1;
}

// 清除所有的覆盖标记
void clear_covers(Munkres* munkres) {
    for (int i = 0; i < munkres->n; i++) {
        munkres->row_covered[i] = false;
        munkres->col_covered[i] = false;
    }
}

// 清除所有的标记零
void erase_primes(Munkres* munkres) {
    for (int i = 0; i < munkres->n; i++) {
        for (int j = 0; j < munkres->n; j++) {
            if (munkres->marked[i][j] == PRIMED) {
                munkres->marked[i][j] = UNMARKED;
            }
        }
    }
}

// 查找矩阵中最小的未覆盖值
float find_smallest(Munkres* munkres) {
    float minval = FLT_MAX;
    float epsilon = 1e-6;
    for (int i = 0; i < munkres->n; i++) {
        if (!munkres->row_covered[i]) {
            for (int j = 0; j < munkres->n; j++) {
                if (!munkres->col_covered[j] && munkres->C[i][j] < minval && munkres->original_C[i][j] != DISALLOWED_VAL) {
                    minval = munkres->C[i][j];
                }
            }
        }
    }
    return minval;
}

// 构建增广路径并调整标记
void convert_path(Munkres* munkres, int count) {
    for (int i = 0; i <= count; i++) {
        int r = munkres->path[i][0];
        int c = munkres->path[i][1];
        if (munkres->marked[r][c] == STARRED) {
            munkres->marked[r][c] = UNMARKED;
        } else {
            munkres->marked[r][c] = STARRED;
        }
    }
}

// 扩展增广路径
void step5_build_path(Munkres* munkres, int* step) {
    int count = 0;
    munkres->path[count][0] = munkres->Z0_r;
    munkres->path[count][1] = munkres->Z0_c;
    count++;

    while (1) {
        int row = find_star_in_col(munkres, munkres->path[count - 1][1]);
        if (row == -1) {
            break;
        }
        munkres->path[count][0] = row;
        munkres->path[count][1] = munkres->path[count - 1][1];
        count++;

        int col = find_prime_in_row(munkres, munkres->path[count - 1][0]);
        if (col == -1) {
            break;
        }
        munkres->path[count][0] = munkres->path[count - 1][0];
        munkres->path[count][1] = col;
        count++;
    }

    // 执行路径转换
    convert_path(munkres, count - 1);

    // 清除覆盖标记和标记零
    clear_covers(munkres);
    erase_primes(munkres);

    *step = 3;
}

// Step 1: 对每一行进行最小值减法
int step1(Munkres* munkres) {
    for (int i = 0; i < munkres->n; i++) {
        float minval = FLT_MAX;
        for (int j = 0; j < munkres->n; j++) {
            if (munkres->C[i][j] < minval && munkres->original_C[i][j] != DISALLOWED_VAL) {
                minval = munkres->C[i][j];
            }
        }
        if (minval == FLT_MAX) {
            // 如果一整行都是DISALLOWED，返回失败状态
            printf("Error: Row %d is entirely DISALLOWED.\n", i);
            return -1;
        }
        for (int j = 0; j < munkres->n; j++) {
            if (munkres->C[i][j] != DISALLOWED_VAL) {
                munkres->C[i][j] -= minval;
            }
        }
    }
    return 2;
}

// Step 2: 标记零
int step2(Munkres* munkres) {
    for (int i = 0; i < munkres->n; i++) {
        for (int j = 0; j < munkres->n; j++) {
            if (fabs(munkres->C[i][j]) < 1e-6 && !munkres->row_covered[i] && !munkres->col_covered[j] && munkres->original_C[i][j] != DISALLOWED_VAL) {
                munkres->marked[i][j] = STARRED;
                munkres->row_covered[i] = true;
                munkres->col_covered[j] = true;
                break;
            }
        }
    }
    clear_covers(munkres);
    return 3;
}

// Step 3: 覆盖包含星号零的所有列
int step3(Munkres* munkres) {
    int count = 0;
    for (int i = 0; i < munkres->n; i++) {
        for (int j = 0; j < munkres->n; j++) {
            if (munkres->marked[i][j] == STARRED && !munkres->col_covered[j]) {
                munkres->col_covered[j] = true;
                count++;
            }
        }
    }

    printf("After Step 3:\n");
    printf("Row covers: ");
    for (int i = 0; i < munkres->n; i++) {
        printf("%d ", munkres->row_covered[i]);
    }
    printf("\nColumn covers: ");
    for (int j = 0; j < munkres->n; j++) {
        printf("%d ", munkres->col_covered[j]);
    }
    printf("\n\n");

    if (count >= munkres->n) {
        return 7; // DONE
    } else {
        return 4;
    }
}

// Step 4: 找到未覆盖的零并标记
int step4(Munkres* munkres) {
    while (find_a_zero(munkres, &munkres->Z0_r, &munkres->Z0_c)) {
        munkres->marked[munkres->Z0_r][munkres->Z0_c] = PRIMED;
        int star_col = find_star_in_row(munkres, munkres->Z0_r);
        if (star_col != -1) {
            munkres->row_covered[munkres->Z0_r] = true;
            munkres->col_covered[star_col] = false;
        } else {
            // 找到一个没有星号零的行，进入Step 5
            return 5;
        }
    }
    // 没有未覆盖的零，进入Step 6
    return 6;
}

// Step 5: 构建增广路径并调整标记
int step5(Munkres* munkres, int* step) {
    step5_build_path(munkres, step);
    return *step;
}

// Step 6: 调整矩阵元素
int step6(Munkres* munkres) {
    float minval = find_smallest(munkres);
    if (minval == FLT_MAX) {
        // 无法调整，矩阵不可解
        printf("Error: Matrix cannot be solved!\n");
        return -1;
    }

    for (int i = 0; i < munkres->n; i++) {
        for (int j = 0; j < munkres->n; j++) {
            if (munkres->original_C[i][j] == DISALLOWED_VAL) {
                continue; // 跳过DISALLOWED位置
            }
            if (munkres->row_covered[i]) {
                munkres->C[i][j] += minval;
            }
            if (!munkres->col_covered[j]) {
                munkres->C[i][j] -= minval;
            }
        }
    }

    // 更新标签
    for (int i = 0; i < munkres->n; i++) {
        if (munkres->row_covered[i]) {
            munkres->lx[i] -= minval;
        }
        for (int j = 0; j < munkres->n; j++) {
            if (munkres->col_covered[j]) {
                munkres->ly[j] += minval;
            }
        }
    }

    printf("After Step 6:\n");
    printf("Row labels (lx): ");
    for (int i = 0; i < munkres->n; i++) {
        printf("%.4lf ", munkres->lx[i]);
    }
    printf("\nColumn labels (ly): ");
    for (int j = 0; j < munkres->n; j++) {
        printf("%.4lf ", munkres->ly[j]);
    }
    printf("\n\n");

    return 4;
}

// 执行Munkres算法并返回状态
int compute(Munkres* munkres) {
    int step = 1;
    while (step != 7) { // 7 是 DONE
        switch (step) {
            case 1:
                step = step1(munkres);
                break;
            case 2:
                step = step2(munkres);
                break;
            case 3:
                step = step3(munkres);
                break;
            case 4:
                step = step4(munkres);
                break;
            case 5:
                step = step5(munkres, &step); // 正确传递步骤变量的地址
                break;
            case 6:
                step = step6(munkres);
                break;
            default:
                printf("Error: Invalid step %d.\n", step);
                return -1;
        }
        if (step == -1) {
            return -1; // 匹配失败
        }
    }
    return 0; // 匹配成功
}

// 获取配对结果
int get_results(Munkres* munkres, Assignment results[], int original_rows, int original_cols) {
    int count = 0;
    for (int i = 0; i < original_rows; i++) {
        for (int j = 0; j < original_cols; j++) {
            if (munkres->marked[i][j] == STARRED && munkres->original_C[i][j] != DISALLOWED_VAL) {
                results[count].row = i;
                results[count].col = j;
                count++;
            }
        }
    }
    return count;
}

// 计算总成本基于原始成本矩阵
float calculate_total_cost(Munkres* munkres, Assignment results[], int count) {
    float total = 0.0;
    for (int i = 0; i < count; i++) {
        int r = results[i].row;
        int c = results[i].col;
        if (munkres->original_C[r][c] != DISALLOWED_VAL) {
            total += munkres->original_C[r][c];
        }
    }
    return total;
}

// 封装的匹配函数
int hungarian_match(float input_matrix[][MAX_SIZE], int input_rows, int input_cols, Assignment results[], int* result_count, float* total_cost) {
    Munkres munkres;
    pad_matrix(&munkres, input_matrix, input_rows, input_cols);
    initialize(&munkres);

    // 打印成本矩阵
    print_matrix(munkres.C, munkres.n, munkres.n, "Cost matrix:");

    // 执行算法
    int status = compute(&munkres);
    if (status != 0) {
        // 匹配失败
        return -1;
    }

    // 获取结果
    *result_count = get_results(&munkres, results, input_rows, input_cols);

    // 计算总成本基于原始成本矩阵
    *total_cost = calculate_total_cost(&munkres, results, *result_count);

    return 0; // 匹配成功
}


// 矩阵取反函数
void invert_matrix(float matrix[MAX_SIZE][MAX_SIZE], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
           if (isinf(matrix[i][j])) {
                matrix[i][j] = 0.0;  // 无穷大（正或负）改为0
            }
            // 如果是正数，变成负数
            else if (matrix[i][j] > 0.0) {
                matrix[i][j] = -matrix[i][j];  // 取反
            }
            // 如果是负数，变成正数
            else if (matrix[i][j] < 0.0) {
                matrix[i][j] = -matrix[i][j];  // 取反
            }
        }
    }
}

// 定义所有测试用例
#define NUM_TESTS 12  // 更新为12个测试用例

int main() {
    // 定义所有测试矩阵和预期结果
    TestCase tests[NUM_TESTS] = {
        // 1. Square
        {
            .matrix = {
                {400, 150, 400},
                {400, 450, 600},
                {300, 225, 300}
            },
            .rows = 3,
            .cols = 3,
            .expected_cost = -1225.0000
        },
        // 2. Rectangular variant
        {
            .matrix = {
                {400, 150, 400, 1},
                {400, 450, 600, 2},
                {300, 225, 300, 3}
            },
            .rows = 3,
            .cols = 4,
            .expected_cost = -1225.0000
        },
        // 3. Square
        {
            .matrix = {
                {10, 10, 8},
                {9, 8, 1},
                {9, 7, 4}
            },
            .rows = 3,
            .cols = 3,
            .expected_cost = -25.0000
        },
        // 4. Square variant with floating point value
        {
            .matrix = {
                {10.1, 10.2, 8.3},
                {9.4, 8.5, 1.6},
                {9.7, 7.8, 4.9}
            },
            .rows = 3,
            .cols = 3,
            .expected_cost = -26.5000
        },
        // 5. Rectangular variant
        {
            .matrix = {
                {10, 10, 8, 11},
                {9, 8, 1, 1},
                {9, 7, 4, 10}
            },
            .rows = 3,
            .cols = 4,
            .expected_cost = -29.0000
        },
        // 6. Rectangular variant with floating point value
        {
            .matrix = {
                {10.01, 10.02, 8.03, 11.04},
                {9.05, 8.06, 1.07, 1.08},
                {9.09, 7.10, 4.11, 10.12}
            },
            .rows = 3,
            .cols = 4,
            .expected_cost = -29.1900
        },
        // 7. Rectangular with DISALLOWED
        {
            .matrix = {
                {4, 5, 6, 0},
                {1, 9, 12, 11},
                {0, 5, 4, 0},
                {12, 12, 12, 10}
            },
            .rows = 4,
            .cols = 4,
            .expected_cost = -34.0000
        },
        // 8. Rectangular variant with DISALLOWED and floating point value
        {
            .matrix = {
                {4.001, 5.002, 6.003, 0},
                {1.004, 9.005, 12.006, 11.007},
                {0, 5.008, 4.009, 0},
                {12.01, 12.011, 12.012, 10.013}
            },
            .rows = 4,
            .cols = 4,
            .expected_cost = -34.0280
        },
        // 9. DISALLOWED to force pairings
        {
            .matrix = {
                {1, 0, 0, 0},
                {0, 2, 0, 0},
                {0, 0, 3, 0},
                {0, 0, 0, 4}
            },
            .rows = 4,
            .cols = 4,
            .expected_cost = -10.0000
        },
        // 10. DISALLOWED to force pairings with floating point value
        {
            .matrix = {
                {1.1, 0, 0, 0},
                {0, 2.2, 0, 0},
                {0, 0, 3.3, 0},
                {0, 0, 0, 4.4}
            },
            .rows = 4,
            .cols = 4,
            .expected_cost = -11.0000
        },
        // 11. Rectangular variant with negative costs
        {
            .matrix = {
                {0.8768, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000},
                {-1.0000,  0.8997, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000},
                {-1.0000, -1.0000,  0.8312, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000},
                {-1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000},
                {-1.0000, -1.0000, -1.0000,  0.8771, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  0.3786,  0.3098, -1.0000,  0.2441, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000},
                {-1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  0.8956,  0.5149, -1.0000, -1.0000, -1.0000, -1.0000,  0.3389, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000},
                {-1.0000, -1.0000, -1.0000, -1.0000,  0.8140, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000}
            },
            .rows = 7,
            .cols = 22,
            .expected_cost = -4.1944
        },
        // 12. Rectangular variant with incomplete columns (last test case)
        {
            .matrix = {
                {0.8768, -1.0000},
                {-1.0000,  0.8997},
                {-1.0000, -1.0000},
                {-1.0000, -1.0000},
                {-1.0000, -1.0000},
                {-1.0000, -1.0000},
                {-1.0000, -1.0000}
            },
            .rows = 7,
            .cols = 2,
            .expected_cost = -1.7765
        }
    };

    // 定义预期的匹配结果（可选）
    // 此处仅验证总成本，匹配结果可根据需要添加

    for (int t = 0; t < NUM_TESTS; t++) {
        TestCase current_test = tests[t];
        printf("=== Test Case %d ===\n", t + 1);

        // 如果需要最大化问题，可以调用 negate_matrix 函数
        invert_matrix(current_test.matrix, current_test.rows, current_test.cols);

        // 执行匹配
        Assignment results[MAX_SIZE];
        int result_count = 0;
        float total_cost = 0.0;
        int status = hungarian_match(current_test.matrix, current_test.rows, current_test.cols, results, &result_count, &total_cost);

        if (status != 0) {
            printf("匹配失败！\n");
            printf("预期的总成本 = %.4lf\n\n", current_test.expected_cost);
            continue;
        }

        // 打印匹配结果
        printf("匹配结果:\n");
        for (int i = 0; i < result_count; i++) {
            int r = results[i].row;
            int c = results[i].col;
            if (current_test.matrix[r][c] != DISALLOWED_VAL) {
                printf("目标 %d 匹配到观测 %d，成本: %.4lf\n", r, c, current_test.matrix[r][c]);
            }
        }
        printf("计算的总成本 = %.4lf\n", total_cost);
        printf("预期的总成本 = %.4lf\n", current_test.expected_cost);

        // 验证结果
        // 使用一个小的误差范围来比较浮点数
        float epsilon = 1e-3;
        if (fabs(total_cost - current_test.expected_cost) < epsilon) {
            printf("测试通过！\n");
        } else {
            printf("测试失败！预期: %.4lf, 得到: %.4lf\n", current_test.expected_cost, total_cost);
        }
        printf("\n");
    }

    return 0;
}
