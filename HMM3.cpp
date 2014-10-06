#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>

struct matrix
{
	int row;
	int col;
	double** m;

	//constructor
	matrix(double** mat, int rows, int cols)
	{
		row = rows;
		col = cols;
		m = mat;
	}

	//constructor 2
	matrix(std::string in)
	{
		std::istringstream iss;
		iss.str(in);
		iss >> row >> col;

        m = initialize(row,col);

		for(int i=0;i<row;++i)
		{
			for(int j=0;j<col;++j)
			iss >> m[i][j];
		}
	}

    //creted the 2d matrix
	static double** initialize(int rows, int cols)
	{
	    double** temp;
	    temp = (double**)calloc(rows , sizeof(double *));
		for(int i=0 ; i< rows ; ++i)
			temp[i] = (double*)calloc(cols , sizeof(double));
        return temp;
	}

	double get(int r, int c)
	{
		return m[r][c];
	}

	void print()
	{
		for (int i=0;i<row;++i)
		{
			for(int j=0;j<col;++j)
				std::cout << m[i][j] <<" ";
			std::cout << std::endl;
		}
	}

	void writeMatrix()
	{
		std::cout << row << " " << col;
		for(int i=0;i<row;++i) {
			for(int j=0; j<col;++j)
				std::cout << " " << m[i][j];
		}
	}

	matrix transpose()
	{
		int i,j;
		double** temp = initialize(col,row);
		for(i=0;i<row;++i)
        {
            for(j=0;j<col;++j)
                temp[j][i] = m[i][j];
        }
        return matrix(temp,col,row);
	}

	matrix getCol(int index)
	{
		int newRows = row;
		int newCols = 1;

		double** temp = initialize(newRows,newCols);

		for(int i=0;i<newRows;++i)
			temp[i][0] = m[i][index];
		return matrix(temp,newRows,newCols);
	}
	
	matrix getRow(int index)
	{
		int newRows = 1;
		int newCols = col;

		double** temp = initialize(newRows,newCols);

		for(int i=0;i<newCols;++i)
			temp[0][i] = m[index][i];
		return matrix(temp,newRows,newCols);
	}

	matrix operator*(matrix H)
	{
		int newRows = row;
		int newCols = H.col;

		if(col != H.row)
			std::cerr << "Fel. matrix 1: " << row << " " << col << ", matrix 2: " << H.row << " " << H.col << std::endl;

		double** temp;
		temp = initialize(newRows,newCols);

		//RÄKNA!
		double t;
		int i,j,k;
		for (i=0;i<newRows;++i)
		{
			for(j=0;j<newCols;++j)
			{
				t = 0;
				//row in first matrix X col in second matrix;
				for(k=0;k<col;++k)
					t+= (m[i][k] * H.get(k,j));
				temp[i][j] = t;
			}
		}
		return matrix(temp,newRows,newCols);
	}

	matrix operator&(matrix H) // element multiplication
	{
		int newRows = row;
		int newCols = col;

		if(row!=H.row || col != H.col)
			std::cout << "FEL " << row << " " << H.row << " , " << col <<" " <<H.col <<std::endl;

		double** temp = initialize(newRows,newCols);

		//RÄKNA!
		int i,j;
		for (i=0;i<newRows;++i)
		{
			for(j=0;j<newCols;++j)
				temp[i][j] = m[i][j]*H.get(i,j);
		}
		return matrix(temp,newRows,newCols);
	}

};

int main(int argc, char **argv)
{
	// Read the file
	std::vector<std::string> board;
	for (std::string line; std::getline(std::cin, line);)
		board.push_back(line);

	matrix A = matrix(board[0]);
	matrix B = matrix(board[1]);
	matrix q = matrix(board[2]);

	std::vector<int> stateSequence;
	int Nstates, index;
	std::istringstream iss;
	iss.str(board[3]);
	iss >> Nstates;
	for(int i=0;i<Nstates;++i) {
		iss >> index;
		stateSequence.push_back(index);
	}

	/**viterbi fast inte samma*/
	
	A.print();
	std::cout << "\n" << std::endl;
	B.print();
	std::cout << "\n" << std::endl;
	
	
	
	int state,i,j,k;
	double maximum = -200;
	//double pr;
	
	double* tempArray = (double*)calloc(A.col,sizeof(double));
	int* tempIndex = (int*)calloc(A.col,sizeof(int));
	
	//xi1	
	int* Sequence = (int*)calloc(Nstates,sizeof(int));
	matrix pr = q&B.getCol(stateSequence[0]).transpose();
	//find max;
	for(i=0;i<pr.col;++i)
	{
		if(pr.get(0,i)> maximum)
		{
			maximum = pr.get(0,i);
			index = i;
		}
	}
	Sequence[0] = index;
	
	for(int s=1;s<Nstates;++s)
	{
		state = stateSequence[s];
		std::cout << "\nstate: " << state << "\nP ";
		pr.print();
		
		for(j=0;j<A.col;++j)
		{
			matrix temp = A.getCol(j).transpose()&pr;
			maximum = -200;
			
			std::cout << "\nTemp: "<< j <<": " << std::endl;
			for(k=0;k<temp.col;++k)
			{
			std::cout << temp.get(0,k)*B.get(j,state) << " (" << k << ") " << std::endl;
				if(temp.get(0,k) > maximum)
				{
					maximum = temp.get(0,k);
					index = k;
				}
			tempArray[j] = maximum*B.get(j,state);
			tempIndex[j] = index;
			}			
		}
		
		//std::cout << "\nMax: "<< tempArray[j] << "\nIndex: " << tempIndex[j] << "\nj = " << j << std::endl;
		
		//find max in tempArray
		maximum = -200;
		for(j=0;j<A.col;++j)
		{
			if(tempArray[j] > maximum)
			{
				maximum = tempArray[j];
				index = tempIndex[j];
				std::cout << "index " << tempIndex[j] << std::endl;
			}
		}
		
		Sequence[s] = index;
		
		//turn tempArrat into a matrix
		double** t =(double**)calloc(1,sizeof(double*));
		t[0] = tempArray;
		pr = matrix(t,1,A.col);
		
	}
	std::cout << "P ";
	pr.print();
	std::cout << "\nSequence ";
	for(i=0;i<Nstates;++i)
		std::cout << Sequence[i] << " ";
	
	
    /**Viterbi algorithm*/
	/*
    int K = A.col;
    int M = B.col;
    int T = Nstates;
    int i,j,k;

    double maximum;
    double temp;

    matrix pi = q&B.getCol(stateSequence[0]).transpose();
    //pi.print();

    double** T1 = A.initialize(K,T);
    double** T2 = A.initialize(K,T);

    for(i=0;i<K;++i)
    {
        T1[i][1] = pi.get(0,i)*B.get(i,stateSequence[0]);
        T2[i][1] = 0;
    }
    for(i=1;i<T;++i)
    {
        for(j=0;j<K;++j)
        {
            maximum = -200;
            index = 0;

            for(k=0;k<K;++k)
            {
                temp = T1[k][i-1]*A.get(k,j)*B.get(j,stateSequence[i]);
                //std::cout << temp << std::endl;
                if(temp>maximum)
                {
                    maximum = temp;
                    index = k;
                }
            }
            T1[j][i] =maximum;
            T2[j][i] =index;
        }
    }
    //double* Z = (double*)calloc(T , sizeof(int));
    int* X = (int*)calloc(T , sizeof(int));

    maximum = -200;
    index = 0;

    for(k=0;k<K;++k)
    {
        if(T1[k][T-1]>maximum)
        {
            maximum = T1[k][T-1];
            index = k;
        }
    }
    //Z[T-1] = index;
    X[T-1] = index;

    for(i=T-1;i>1;--i)
    {
        //Z[i-1] = T2[Z[i][i];
        X[i-1] = T2[X[i]][i];
    }

    for(i=T-1;i>=0;--i)
        std::cout << X[i]<< " ";
    std::cout<<std::endl;
	*/

	return 0;
}

