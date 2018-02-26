/*
Primal���⣺ ����xT = [x1,x2,x3,x4,x5,x6]
min: -7x1 + 7x2 - 2x3 - x4 - 6x5
s.t.:
    3x1 - x2 + x3 - 2x4 = -3
    2x1 + x2 + x4 + x5 = 4
    -x1 + 3x2 - 3x4 + x6 = 12
    xi >= 0; (i = 1,...,6)

Dual���⣺  ����yT = [y1,y2,y3]
max:-3y1 + 4y2 + 12y3
s.t.:
    3y1 + 2y2 - y3 <= -7
    -y1 + y2 + 3y3 <= 7
    y1 + y2 <= -2
    -2y1 -2y2 -3y3 <= -1
    y2 <= -6
    y3 <= 0

��ż���������ļ���7.txt
*/

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
   double ** matrix;//LP���� ����ԭ����Ŀ�꺯��ϵ��c  s.t.ϵ��A�ͽ��b
   char * sign;// s.t.���� ���ڵ��ڣ�-1  С�ڵ��ڣ�1  ���ڣ�0
   char stop;//solution״̬ 1���п��н� 2����ΪUnbounded 3���޿��н�
   int * state;//matrix����ֵ���ڻ�������stateΪ0 �ڻ�������Ϊ�������±�i
   int * base;//��ǰ��
   int ops;//��������
   int col;//������-��ż��������
   int row;//Լ����-��ż��������
   int n1;
   int m1;
}task;


void readFile(task * current, char * filename)//��ȡLP�����ļ�
{
        int i, j, no, rr;
	    double ch;

        FILE * fdata = fopen(filename, "r");

        fflush(stdin);
        fscanf(fdata, "%d %d\n", &(current->col), &(current->row));   //n��m��  n=6 m=3 ԭ����6����3Լ��
        current->n1 = current->col + 1;
        current->m1 = current->row + 1;
        current->matrix = (double**)malloc((current->row + 2) * sizeof(double*));
        current->sign = (char*)malloc((current->row + 1)*sizeof(char));
        if ( current->matrix == NULL || current->sign == NULL)
        		exit(-1);

        //malloc+��ʼ��matrix
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

        //��ʼ����LP�ļ�����
		for(i = 1; i <= current->row; i++)
			fscanf(fdata, "%lf", &(current->matrix[i][0]));  //��ȡm������  -3.00	4.00	12.00

		for(i = 1; i <= current->row; i++)
			fscanf(fdata, "%hd", &(current->sign[i])); //��ȡm������0	0	0  ��P����AX=b����������

		// ��ȡϵ������
		for (i = 1; i <= current->col; i++) //n��   n=6
		{
			//ʾ����-7	3	1	3	2	2	3	-1  ��3y1 + 2y2 -y3 <= -7
            fscanf (fdata, "%lf", &ch); //-7
            current->matrix[0][i] = -ch;//7

			// ��ȡ��0ϵ��
            fscanf (fdata, "%d", &no);//3��
            for (j = 1; j <= no; j++)
            {
                fscanf (fdata, "%d", &rr );//1 2 3
                fscanf (fdata, "%lf", &(current->matrix[rr][i])); //3  2  -1
            }

		}
        fclose(fdata);
}

void insertSlacks(task * current)//ԭ���������ɳڱ���
{
    int i, j;
    for (i = 1; i <= current->row; i++)//�ɳڱ�����=Լ����
    {
        if (current->sign[i] == 0) continue;// =
        current->base[i] = current->n1;//������ɳڱ�������Ȼ��һ���->����base
        current->state[current->n1] = i;//����state
        current->matrix[i][current->n1] = (double)current->sign[i];//����matrix
        current->n1++;//���ӱ�����
        if (current->sign[i] == -1) // >=
            for (j = 0; j < current->n1; j++)
                current->matrix[i][j] *= -1;
    }
}

void findbase(task * current) //��ʼ��
{
    int i, j;
    for (j = current->n1 - 1; j > 0; j--)
    {
        if (current->state[j] != 0) continue;//state�Ѹ���
        for (i = current->m1 - 1; i > 0; i--)
        {
            if (current->base[i] != 0) continue;//�Ѿ��ǻ�����
            if (current->matrix[i][j] != 0)//�ǻ� matrix[i][j]��0
            {
                Pivot(current,i,j);//�б任
                current->base[i] = j; //����������base
                current->state[j] = i;//���и���state
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

    current->stop = 1;//matrix[i][0]��>=0����biû�и���->�ҵ�solution��ֹͣ
    return 1;
}

int in_base(task * current)  //�����
{
    int i, row = 1;
    for (i = 2; i < current->m1; i++)
        if (current->matrix[i][0] < current->matrix[row][0])//δstop->bi�и��� ѡbi��С��i���
            row = i;//��¼���������Ӧ��row �����ҳ���

    return row;
}

int out_base(task * current, int row) //�ҳ���
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
        if (current->matrix[0][j] / current->matrix[row][j] <= rMin)// cj/aij��С��j��Ӧ������->����e
        {
            rMin = current->matrix[0][j] / current->matrix[row][j];
            col = j;
        }
    }

    if (count == current->n1 - 1)
    {
        current->stop = 3;//��e=���� �޿��н�
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
            current->matrix[row][j] /= pivot; //�������������->1

    //�µĻ��䵥λ��
    for (i = 0; i < current->m1; i++)
    {
        if (i == row) continue;
		if (current->matrix[i][col] == 0) continue;
        t = current->matrix[i][col];
        for (j = 0; j < current->n1; j++)
            current->matrix[i][j] -= current->matrix[row][j] * t;//�б任
    }
}

void print(task * current)//���solution
{
    int j;
    FILE *fresults = fopen ("result.txt", "w");//����ļ�

    fflush(fresults);
    if (current->stop >= 2)
    {
        if (current->stop == 2)//solutionΪ����
            fprintf ( fresults, "Cost: Unbounded\n");
        if (current->stop == 3)//no solution
            fprintf ( fresults, "Cost: Infeasible\n");
    }
    else
    {
        fprintf (fresults, "Cost: %.2lf \n", current->matrix[0][0]);//Ŀ�꺯��ֵ
        fprintf (fresults, "solution:\n");
        for(j = 1; j < current->n1; j++)//solution
            if (current->state[j] != 0)
                fprintf(fresults, "x%d  = %.2lf\n", j, current->matrix[current->state[j]][0]);
            else
                fprintf(fresults, "x%d  = %.2lf\n", j, 0);
    }
    fclose(fresults);
}


void dualSimplex(task * current)//�㷨����
{
        int i, row, col;

		//��ʼ��
        current->ops = 0;
        current->stop = 0;
        current->state = (int*)malloc((current->col + current->row + 2) * sizeof(int));
        current->base = (int*)malloc((current->row + 2)*sizeof(int));

		for(i = 0; i < current->col + current->row + 2; i++)
   			current->state[i] = 0;

		for(i = 0;i < current->row + 2; i++)
   			current->base[i] = 0;

        //ԭ���������ɳڱ���
        //insertSlacks(current);
        //�����ʼ��
        findbase(current);

        //ԭ���ⲻ����  һֱimprove
        while (!stopping(current))  //biȫ����0 stop
        {
            current->ops++;
            row = in_base(current);//�����

            //�ҳ���
            col = out_base(current,row);
            if (col == 0) return;  //Infeasible

            //���»�->base state
            current->state[current->base[row]] = 0;
            current->state[col] = row;
            current->base[row] = col;

            //���»�->��ԭLP������˹�б任
            Pivot(current,row,col);

         }
}

void Free (task *current)//free��matrix���� state��base��state����
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
