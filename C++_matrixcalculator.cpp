#include <iostream>
#include <conio.h>
using namespace std;

class MatrixException
{
private:
    std::string m_error;

public:
    MatrixException(std::string error)
        : m_error(error)
    {
    }

    const char* getError()
    {
        return m_error.c_str();
    }
};

class nonSquareMatrix
{
private:
    unsigned int numberRows, numberColumns;
    double **matrix;
public:
    nonSquareMatrix &operator=(const nonSquareMatrix &mat);
    void swapRows(int row1, int row2);
    void multiplyRow(int row, double factor);
    nonSquareMatrix reducedRowEchelonForm() const;
    void setElement(int row, int column, double value);
    double getElement(int row, int column);
    void setDimensions(int rows, int columns);
    friend ostream &operator << (ostream &cout, const nonSquareMatrix &mat);
    friend istream &operator >> (istream &cin, nonSquareMatrix &mat);
};
nonSquareMatrix &nonSquareMatrix::operator=(const nonSquareMatrix &mat)
{
    if(numberRows!=0)
    {
        for (int i=0; i<numberRows; ++i)
            delete [] matrix[i];
        delete [] matrix;
    }
    numberColumns=mat.numberColumns;
    numberRows=mat.numberRows;
    matrix=new double*[numberRows];
    for(int i=0; i<numberRows; ++i)
        matrix[i]=new double[numberColumns];
    for(int i=0; i<numberRows; ++i)
        for(int j=0; j<numberColumns; ++j)
            matrix[i][j]=mat.matrix[i][j];
    return *this;
}

void nonSquareMatrix::swapRows(int row1, int row2)
{
    for(int i=0; i<numberColumns; ++i)
        swap(matrix[row1][i],matrix[row2][i]);
}

void nonSquareMatrix::multiplyRow(int row, double factor)
{
    for(int i=0; i<numberColumns; ++i)
        matrix[row][i]*=factor;
}

nonSquareMatrix nonSquareMatrix::reducedRowEchelonForm() const
{
    nonSquareMatrix mat=*this;
    for(int i=0; i<min(mat.numberColumns, mat.numberRows); ++i)
    {
        int firstNonZeroRow=-1;
        for(int j=i; j<mat.numberRows; ++j)
            if(mat.matrix[j][i]!=0)
            {
                firstNonZeroRow=j;
                break;
            }
        if(firstNonZeroRow!=i)
            mat.swapRows(i, firstNonZeroRow);
        mat.multiplyRow(i,(1/mat.matrix[i][i]));
        for(int j=i+1; j<mat.numberRows; ++j)
        {
            double coefficient=(-1)*mat.matrix[j][i];
            for(int q=i; q<mat.numberColumns; ++q)
                mat.matrix[j][q]+=coefficient*mat.matrix[i][q];
        }
    }
    for(int i=mat.numberRows-1; i>=0; --i)
        if(mat.matrix[i][i]==1)
            for(int j=i-1; j>=0; --j)
            {
                double coefficient=(-1)*mat.matrix[j][i];
                for(int q=i; q<mat.numberColumns; ++q)
                    mat.matrix[j][q]+=coefficient*mat.matrix[i][q];
            }
    return mat;
}

void nonSquareMatrix::setElement(int row, int column, double value)
{
    matrix[row][column]=value;
}

void nonSquareMatrix::setDimensions(int rows, int columns)
{
    numberColumns=columns, numberRows=rows;
    matrix=new double*[numberRows];
    for(int i=0; i<numberRows; ++i)
        matrix[i]=new double[numberColumns];
    for(int i=0; i<numberRows; ++i)
        for(int j=0; j<numberColumns; ++j)
            matrix[i][j]=0;
}

double nonSquareMatrix::getElement(int row, int column)
{
    return matrix[row][column];
}
istream &operator >> (istream &cin, nonSquareMatrix &mat)
{
    cin>>mat.numberRows>>mat.numberColumns;
    mat.matrix=new double*[mat.numberRows];
    for(int i=0; i<mat.numberRows; ++i)
        mat.matrix[i]=new double[mat.numberColumns];
    for(int i=0; i<mat.numberRows; ++i)
        for(int j=0; j<mat.numberColumns; ++j)
            cin>>mat.matrix[i][j];
    return cin;
}
ostream &operator << (ostream &cout, const nonSquareMatrix &mat)
{
    for(int i=0; i<mat.numberRows; ++i)
    {
        for(int j=0; j<mat.numberColumns; ++j)
            cout<<mat.matrix[i][j]<<" ";
        cout<<"\n";
    }
    cout<<"\n";
    return cout;
}

class Matrix
{
private:
    unsigned int order;
    double **matrix;
public:
    Matrix();
    Matrix(int newOrder, double **newMatrix);
    Matrix(const Matrix &mat);
    ~Matrix();
    Matrix &operator=(const Matrix &mat);
    Matrix operator+(const Matrix &mat);
    Matrix operator-(const Matrix &mat);
    Matrix operator*(const Matrix &mat);
    bool operator==(const Matrix &mat);
    bool operator!=(const Matrix &mat);
    Matrix operator&(const double &scalar);
    Matrix operator^(const int &power);
    double determinant() const;
    void swapRows(int row1, int row2);
    void setOrder(int newOrder);
    Matrix inverseMatrix() const;
    friend ostream &operator << (ostream &cout, const Matrix &mat);
    friend istream &operator >> (istream &cin, Matrix &mat);
};
Matrix::Matrix()
{
    order=0;
}
Matrix::Matrix(int newOrder, double **newMatrix)
{
    order=newOrder;
    matrix=new double*[order];
    for(int i=0; i<order; ++i)
        matrix[i]=new double[order];
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            matrix[i][j]=newMatrix[i][j];
}
Matrix::Matrix(const Matrix &mat)
{
    order=mat.order;
    matrix=new double*[order];
    for(int i=0; i<order; ++i)
        matrix[i]=new double[order];
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            matrix[i][j]=mat.matrix[i][j];
}
Matrix::~Matrix()
{
    for (int i=0; i<order; ++i)
        delete [] matrix[i];
    delete [] matrix;
    order=0;
}
Matrix &Matrix::operator=(const Matrix &mat)
{
    if(order!=0)
    {
        for (int i=0; i<order; ++i)
            delete [] matrix[i];
        delete [] matrix;
    }
    order=mat.order;
    matrix=new double*[order];
    for(int i=0; i<order; ++i)
        matrix[i]=new double[order];
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            matrix[i][j]=mat.matrix[i][j];
    return *this;
}
Matrix Matrix::operator+(const Matrix &mat)
{
    if(order!=mat.order)
        throw MatrixException("Matrices are of different orders");
    Matrix result=*this;
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            result.matrix[i][j]+=mat.matrix[i][j];
    return result;
}
Matrix Matrix::operator-(const Matrix &mat)
{
    if(order!=mat.order)
        throw MatrixException("Matrices are of different orders");
    Matrix result=*this;
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            result.matrix[i][j]-=mat.matrix[i][j];
    return result;
}
Matrix Matrix::operator*(const Matrix &mat)
{
    if(order!=mat.order)
        throw MatrixException("Matrices are of different orders");
    Matrix result;
    result.setOrder(order);
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            for(int k=0; k<order; ++k)
                result.matrix[i][j]+=matrix[i][k]*mat.matrix[k][j];
    return result;
}
bool Matrix::operator==(const Matrix &mat)
{
    if(order!=mat.order)
        return false;
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            if(matrix[i][j]!=mat.matrix[i][j])
                return false;
    return true;
}
bool Matrix::operator!=(const Matrix &mat)
{
    return (!(*this==mat));
}
Matrix Matrix::operator&(const double &scalar)
{
    Matrix result=*this;
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            result.matrix[i][j]*=scalar;
    return result;
}
Matrix Matrix::operator^(const int &power)
{
    Matrix result;
    if(power==0)
    {
        result.setOrder(order);
        for(int i=0; i<order; ++i)
            result.matrix[i][i]=1;
        return result;
    }
    if(power<0)
    {
        Matrix inverse=inverseMatrix();
        result=inverse;
        int pow=power*(-1);
        for(int i=2; i<=pow; ++i)
            result=result*inverse;
        return result;
    }
    result=*this;
    for(int i=2; i<=power; ++i)
        result=result*(*this);
    return result;
}
double Matrix::determinant() const
{
    Matrix mat=*this;
    double matrixDeterminant=1;
    for(int i=0; i<mat.order-1; ++i)
    {
        int firstNonZeroRow=-1;
        for(int j=i; j<mat.order; ++j)
            if(mat.matrix[j][i]!=0)
            {
                firstNonZeroRow=j;
                break;
            }
        if(firstNonZeroRow!=i)
            matrixDeterminant*=(-1), mat.swapRows(i, firstNonZeroRow);
        for(int j=i+1; j<order; ++j)
        {
            double coefficient=(-1)*mat.matrix[j][i]/mat.matrix[i][i];
            for(int q=i; q<order; ++q)
                mat.matrix[j][q]+=coefficient*mat.matrix[i][q];
        }
    }
    for(int i=0; i<order; ++i)
        matrixDeterminant*=mat.matrix[i][i];
    return matrixDeterminant;
}
void Matrix::setOrder(int newOrder)
{
    if(order!=0)
    {
        for (int i=0; i<order; ++i)
            delete [] matrix[i];
        delete [] matrix;
    }
    order=newOrder;
    matrix=new double*[order];
    for(int i=0; i<order; ++i)
        matrix[i]=new double[order];
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            matrix[i][j]=0;
}
void Matrix::swapRows(int row1, int row2)
{
    for(int i=0; i<order; ++i)
        swap(matrix[row1][i],matrix[row2][i]);
}
Matrix Matrix::inverseMatrix() const
{
    if(determinant()==0)
        throw MatrixException("Matrix is not invertible");
    Matrix inverse=*this;
    nonSquareMatrix auxInverse;
    auxInverse.setDimensions(order, 2*order);
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            auxInverse.setElement(i,j, matrix[i][j]);
    for(int i=0; i<order; ++i)
        auxInverse.setElement(i,order+i,1);
    auxInverse=auxInverse.reducedRowEchelonForm();
    for(int i=0; i<order; ++i)
        for(int j=0; j<order; ++j)
            inverse.matrix[i][j]=auxInverse.getElement(i,j+order);
    return inverse;
}
istream &operator >> (istream &cin, Matrix &mat)
{
    cin>>mat.order;
    mat.matrix=new double*[mat.order];
    for(int i=0; i<mat.order; ++i)
        mat.matrix[i]=new double[mat.order];
    for(int i=0; i<mat.order; ++i)
        for(int j=0; j<mat.order; ++j)
            cin>>mat.matrix[i][j];
    return cin;
}
ostream &operator << (ostream &cout, const Matrix &mat)
{
    for(int i=0; i<mat.order; ++i)
    {
        for(int j=0; j<mat.order; ++j)
            cout<<mat.matrix[i][j]<<" ";
        cout<<"\n";
    }
    cout<<"\n";
    return cout;
}
void printOptions1()
{
    cout<<"Select the operation you wish to perform\n";
    cout<<"0.The sum of the two matrices\n";
    cout<<"1.The difference of the two matrices\n";
    cout<<"2.The product of the two matrices\n";
    cout<<"3.Check whether the two matrices are the same\n";
    cout<<"4.Multiply one of the matrices with a scalar\n";
    cout<<"5.Raise one matrix to a certain power\n";
    cout<<"6.The determinant of both matrices\n";
    cout<<"7.The inverse of both matrices\n";
    cout<<"8.Enter a new pair of matrices\n";
    cout<<"9.Exit the program\n";
}
void printOptions2()
{
    cout<<"Select the operation you wish to perform\n";
    cout<<"1.Go back to main menu\n";
    cout<<"2.Exit the program\n";
}
void userMenu()
{
    Matrix matrix1, matrix2;
    cout<<"Enter a positive number n and a matrix of that order:\n";
    cin>>matrix1;
    cout<<"Enter the order(preferably the same as that of the first one)";
    cout<<"and the elements of a second matrix:\n";
    cin>>matrix2;
    while(1)
    {
        system("cls");
        cout<<"The two matrices are:\n";
        cout<<matrix1<<matrix2;
        printOptions1();
        char option;
        option=getch();
        if(option=='0')
        {
            system("cls");
            try
            {
                cout<<"The sum of the two matrices is\n"<<matrix1+matrix2;
            }
            catch (MatrixException &exception)
            {
                std::cerr << "A matrix exception occurred (" << exception.getError() << ")\n";
            }
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='1')
        {
            system("cls");
            try
            {
                cout<<"The difference of the two matrices is\n"<<matrix1-matrix2;
            }
            catch (MatrixException &exception)
            {
                std::cerr << "A matrix exception occurred (" << exception.getError() << ")\n";
            }
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='2')
        {
            system("cls");
            try
            {
                cout<<"The product of the two matrices is\n"<<matrix1*matrix2;
            }
            catch (MatrixException &exception)
            {
                std::cerr << "A matrix exception occurred (" << exception.getError() << ")\n";
            }
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='3')
        {
            system("cls");
            if(matrix1==matrix2)
                cout<<"The two matrices are the same\n";
            if(matrix1!=matrix2)
                cout<<"The two matrices are different\n";
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='4')
        {
            system("cls");
            cout<<"Specify which matrix you would like to multiply\n";
            int orderOfMatrix;
            double scalar;
            cin>>orderOfMatrix;
            cout<<"Enter the scalar:\n";
            cin>>scalar;
            cout<<"The result is:\n";
            if(orderOfMatrix==1)
                cout<<(matrix1&scalar);
            else
                cout<<(matrix2&scalar);
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='5')
        {
            system("cls");
            cout<<"Specify which matrix you would like to raise to a power\n";
            int orderOfMatrix;
            int power;
            cin>>orderOfMatrix;
            cout<<"Enter the power:\n";
            cin>>power;
            cout<<"The result is:\n";
            if(orderOfMatrix==1)
                cout<<(matrix1^power);
            else
                cout<<(matrix2^power);
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='6')
        {
            system("cls");
            cout<<"The determinant of the first matrix is "<<matrix1.determinant()<<"\n";
            cout<<"The determinant of the second matrix is "<<matrix2.determinant()<<"\n";
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='7')
        {
            system("cls");
            try
            {
                cout<<"The inverse of the first matrix is\n"<<matrix1.inverseMatrix();
            }
            catch (MatrixException &exception)
            {
                std::cerr << "A matrix exception occurred for the first matrix(" << exception.getError() << ")\n";
            }
            try
            {
                cout<<"The inverse of the second matrix is\n"<<matrix2.inverseMatrix();
            }
            catch (MatrixException &exception)
            {
                std::cerr << "A matrix exception occurred for the second matrix(" << exception.getError() << ")\n";
            }
            printOptions2();
            char option2;
            option2=getch();
            if(option2=='1')
                continue;
            if(option2=='2')
                exit(0);
            cout<<"Selected option is invalid.\n";
            exit(0);
        }
        if(option=='8')
        {
            system("cls");
            cout<<"Enter the new values for the first matrix:\n";
            cin>>matrix1;
            cout<<"Enter the new values for the second matrix\n";
            cin>>matrix2;
            continue;
        }
        if(option=='9')
            exit(0);
        cout<<"Selected option is invalid.\n";
        exit(0);
    }
}
int main()
{
    userMenu();
    return 0;
}
