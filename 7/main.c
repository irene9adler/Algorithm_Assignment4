/*
Primal问题： 变量xT = [x1,x2,x3,x4,x5,x6]
min: -7x1 + 7x2 - 2x3 - x4 - 6x5
s.t.:
    3x1 - x2 + x3 - 2x4 = -3
    2x1 + x2 + x4 + x5 = 4
    -x1 + 3x2 - 3x4 + x6 = 12
    xi >= 0; (i = 1,...,6)

Dual问题：  变量yT = [y1,y2,y3]
max:-3y1 + 4y2 + 12y3
s.t.:
    3y1 + 2y2 - y3 <= -7
    -y1 + y2 + 3y3 <= 7
    y1 + y2 <= -2
    -2y1 -2y2 -3y3 <= -1
    y2 <= -6
    y3 <= 0

对偶问题描述文件：7.txt
*/

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
   double ** matrix;//LP求解表 包括原问题目标函数系数c  s.t.系数A和结果b
   char * sign;// s.t.符号 大于等于：-1  小于等于：1  等于：0
   char stop;//solution状态 1：有可行解 2：解为Unbounded 3：无可行解
   int * state;//matrix中数值不在基向量中state为0 在基向量中为基向量下标i
   int * base;//当前基
   int ops;//迭代次数
   int col;//变量数-对偶问题列数
   int row;//约束数-对偶问题行数
   int n1;
   int m1;
}task;


void readFile(task * current, char * filename)//读取LP描述文件
{
        int i, j, no, rr;
	    double ch;

        FILE * fdata = fopen(filename, "r");

        fflush(stdin);
        fscanf(fdata, "%d %d\n", &(current->col), &(current->row));   //n行m列  n=6 m=3 原问题6变量3约束
        current->n1 = current->col + 1;
        current->m1 = current->row + 1;
        current->matrix = (double**)malloc((current->row + 2) * sizeof(double*));
        current->sign = (char*)malloc((current->row + 1)*sizeof(char));
        if ( current->matrix == NULL || current->sign == NULL)
        		exit(-1);

        //malloc+初始化matrix
        for (i = 0; i < current->row + 2; i++)
        {
            current->matrix[i] = (double*)malloc((current->col + current->row + 2)*sizeof(double));
            if (current->matrix[i] == NULL)
                exit(-1);
        }

        for (i = 0; i < current->row + 2; i ++)
        {
            for ( j = 0; j < current->col + current->row + 2; j++)
                current->matrix[i][j] = 0.0;
        }

        //开始读入LP文件数据
		for(i = 1; i <= current->row; i++)
			fscanf(fdata, "%lf", &(current->matrix[i][0]));  //读取m列数据  -3.00	4.00	12.00

		for(i = 1; i <= current->row; i++)
			fscanf(fdata, "%hd", &(current->sign[i])); //读取m列数据0	0	0  即P问题AX=b的三个等于

		// 读取系数矩阵
		for (i = 1; i <= current->col; i++) //n行   n=6
		{
			//示例：-7	3	1	3	2	2	3	-1  即3y1 + 2y2 -y3 <= -7
            fscanf (fdata, "%lf", &ch); //-7
            current->matrix[0][i] = -ch;//7

			// 读取非0系数
            fscanf (fdata, "%d", &no);//3项
            for (j = 1; j <= no; j++)
            {
                fscanf (fdata, "%d", &rr );//1 2 3
                fscanf (fdata, "%lf", &(current->matrix[rr][i])); //3  2  -1
            }

		}
        fclose(fdata);
}

void insertSlacks(task * current)//原问题引入松弛变量
{
    int i, j;
    for (i = 1; i <= current->row; i++)//松弛变量数=约束数
    {
        if (current->sign[i] == 0) continue;// =
        current->base[i] = current->n1;//引入的松弛变量是天然的一组基->插入base
        current->state[current->n1] = i;//更新state
        current->matrix[i][current->n1] = (double)current->sign[i];//更新matrix
        current->n1++;//增加变量数
        if (current->sign[i] == -1) // >=
            for (j = 0; j < current->n1; j++)
                current->matrix[i][j] *= -1;
    }
}

void findbase(task * current) //初始基
{
    int i, j;
    for (j = current->n1 - 1; j > 0; j--)
    {
        if (current->state[j] != 0) continue;//state已更新
        for (i = current->m1 - 1; i > 0; i--)
        {
            if (current->base[i] != 0) continue;//已经是基向量
            if (current->matrix[i][j] != 0)//非基 matrix[i][j]非0
            {
                Pivot(current,i,j);//行变换
                current->base[i] = j; //列向量放入base
                current->state[j] = i;//此列更新state
            }
        }
    }
}

char stopping(task * current)
{
    int i;
    for (i = 1; i < current->m1; i++)
        if (current->matrix[i][0] < 0)
            return 0;

    current->stop = 1;//matrix[i][0]均>=0，即bi没有负数->找到solution，停止
    return 1;
}

int in_base(task * current)  //找入基
{
    int i, row = 1;
    for (i = 2; i < current->m1; i++)
        if (current->matrix[i][0] < current->matrix[row][0])//未stop->bi有负数 选bi最小的i入基
            row = i;//记录入基向量对应的row 用于找出基

    return row;
}

int out_base(task * current, int row) //找出基
{
    int j, col = current->n1 - 1;
    int count = 0;
    double rMin = 1000000;

    for (j = current->n1 - 1; j > 0; j--)
    {
        if (current->matrix[row][j] > 0)//
        {
            count++;
            continue;
        }
        if (current->matrix[0][j] / current->matrix[row][j] <= rMin)// cj/aij最小的j对应列向量->出基e
        {
            rMin = current->matrix[0][j] / current->matrix[row][j];
            col = j;
        }
    }

    if (count == current->n1 - 1)
    {
        current->stop = 3;//△e=无穷 无可行解
        return 0;
    }
    return col;
}

void Pivot(task * current,int row,int col)
{
    int i, j;
    double t, pivot = current->matrix[row][col];

    if (pivot != 1)
        for (j = 0; j < current->n1; j++)
            current->matrix[row][j] /= pivot; //新入基的列向量->1

    //新的基变单位阵
    for (i = 0; i < current->m1; i++)
    {
        if (i == row) continue;
		if (current->matrix[i][col] == 0) continue;
        t = current->matrix[i][col];
        for (j = 0; j < current->n1; j++)
            current->matrix[i][j] -= current->matrix[row][j] * t;//行变换
    }
}

void print(task * current)//输出solution
{
    int j;
    FILE *fresults = fopen ("result.txt", "w");//结果文件

    fflush(fresults);
    if (current->stop >= 2)
    {
        if (current->stop == 2)//solution为无穷
            fprintf ( fresults, "Cost: Unbounded\n");
        if (current->stop == 3)//no solution
            fprintf ( fresults, "Cost: Infeasible\n");
    }
    else
    {
        fprintf (fresults, "Cost: %.2lf \n", current->matrix[0][0]);//目标函数值
        fprintf (fresults, "solution:\n");
        for(j = 1; j < current->n1; j++)//solution
            if (current->state[j] != 0)
                fprintf(fresults, "x%d  = %.2lf\n", j, current->matrix[current->state[j]][0]);
            else
                fprintf(fresults, "x%d  = %.2lf\n", j, 0);
    }
    fclose(fresults);
}


void dualSimplex(task * current)//算法主体
{
        int i, row, col;

		//初始化
        current->ops = 0;
        current->stop = 0;
        current->state = (int*)malloc((current->col + current->row + 2) * sizeof(int));
        current->base = (int*)malloc((current->row + 2)*sizeof(int));

		for(i = 0; i < current->col + current->row + 2; i++)
   			current->state[i] = 0;

		for(i = 0;i < current->row + 2; i++)
   			current->base[i] = 0;

        //原问题引入松弛变量
        //insertSlacks(current);
        //计算初始基
        findbase(current);

        //原问题不可行  一直improve
        while (!stopping(current))  //bi全大于0 stop
        {
            current->ops++;
            row = in_base(current);//找入基

            //找出基
            col = out_base(current,row);
            if (col == 0) return;  //Infeasible

            //更新基->base state
            current->state[current->base[row]] = 0;
            current->state[col] = row;
            current->base[row] = col;

            //更新基->对原LP表做高斯行变换
            Pivot(current,row,col);

         }
}

void Free (task *current)//free掉matrix矩阵 state、base、state数组
{
    free(current->state);
    free(current->base);
    free(current->sign);

    int i;
    for (i = 0; i < current->row + 2; i++)
        free (current->matrix[i]);
    free(current->matrix);
}


int main()
{
    task Task;
    readFile(&Task,"7.txt");
    dualSimplex(&Task);
    print(&Task);
    Free(&Task);

    return 0;
}
