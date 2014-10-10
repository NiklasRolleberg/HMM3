#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <limits>

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
	int T,N,M;
	N = A.col;
	M = B.col;
	std::istringstream iss;
	iss.str(board[3]);
	iss >> T;
	for(int i=0;i<T;++i) {
		int index;
		iss >> index;
		stateSequence.push_back(index);
	}

	/**viterbi fast inte samma*/
	/*
	A.print();
	std::cout << std::endl;
	B.print();
	std::cout << std::endl;
	q.print();
	std::cout << "------------\n" << std::endl;
	*/

	int state;//,i,k;
	double maximum = -std::numeric_limits<double>::max();

	double* tempArray = (double*)calloc(N,sizeof(double));
	double* prob = (double*)calloc(N,sizeof(double));
	int sequences[T][N];
    
	//t=0;
	for(int i=0;i<N;++i)
	{
		prob[i] = q.get(0,i)*B.get(i,stateSequence[0]);
		int index = 0;
		if(prob[i]> maximum)
		{
			maximum = prob[i];
			index = i;
		}
		sequences[0][i] = index;
	}

	/*
	std::cout << "\nstage 0" << std::endl;
	for(int i=0;i<A.col;++i)
		std::cout << prob[i] << " ";
	std::cout<<"\n\n"<<std::endl;
	*/

	for(int t=1;t<T;++t)
	{
		state = stateSequence[t];

		for(int i=0;i<N;++i)
		{
			maximum = prob[0] * A.get(0,i) * B.get(i,state);;
			int index = 0;
			//std::cout << prob[0] << " x "  << A.get(0,i) <<  " x " << B.get(i,state) <<" "<<"("<< 0 << ")"<<"| ";
			for(int k=1;k<N;++k)
			{
				double te = prob[k] * A.get(k,i) * B.get(i,state);
				//std::cout << prob[k] << " x "  << A.get(k,i) <<  " x " << B.get(i,state) <<" "<<"("<< k << ")"<<"| ";
				if(te > maximum)
				{
					maximum = te;
					index = k;
				}
			}
			//std::cout << "\nt: " << t <<"\ni: " << i << "\nindex: " << index << std::endl;
			tempArray[i] = maximum;
			sequences[t][i] = index;
		}

		for(int i=0;i<N;++i)
			prob[i] = tempArray[i];

		/*
		std::cout << "stage "<< t << std::endl;
		for(int i=0;i<B.col;++i)
			std::cout << prob[i] << " ";
		std::cout<<"\n\n"<<std::endl;
		*/
	}
	/*
	for(int t=0;t<T;++t)
	{
		for(int i=0;i<T;++i)
			std::cout << sequences[t][i] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	
	int current = 0;
	maximum = prob[0];
	for(int i=1;i<N;++i)
	{
		if(prob[i]>maximum)
		{
			maximum = prob[i];
			current = i;
		}
	}

	std::vector<int> e;
	//std::cout <<"i " << i <<" current: " << current<< std::endl;
	for(int t = (T-1);t>=0;t--)
	{
		e.push_back(current);
		current = sequences[t][current];
		//std::cout <<"i " << i <<" current: " << current<< std::endl;
	}
	//e.push_back(current);

    //std::cout << e.size() << std::endl;

	for(int i=e.size()-1;i>=0;i--)
    {
        //std::cout << " ("<<i<<")->";
        std::cout << e[i]<< " ";
    }
	std::cout<<std::endl;
	
	
	return 0;
}

