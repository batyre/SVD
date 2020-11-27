#include <iostream>
#include <cmath>
#include <fstream>
#include <istream>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <vector>
//#include <chrono>

using namespace std;

const double PI = 3.141592653589793;

double getRandomNumber(int Min, int Max)
{
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
    // Равномерно распределяем рандомное число в нашем диапазоне
    return static_cast<double>(rand() * fraction * (Max - Min + 1) + Min);
}
//Класс для замера времени выполнения
//class SimpleTimer
//{
//private:
//    chrono::time_point<chrono::steady_clock> start;
//    chrono::time_point<chrono::steady_clock> end;
//public:
//  SimpleTimer()
    //{
      //start= chrono::high_resolution_clock::now();
//  }
//~SimpleTimer()
    //{
      //end= chrono::high_resolution_clock::now();
        //chrono::duration<float> duration = end - start;
        //cout << "Duration = " << duration.count() << " s" << endl;
//  }


//};
// шаблонный класс Матрица
template <typename T>
class MATRIX
{
private:
    T** M; // матрица
    int m; // количество строк
    int n; // количество столбцов

public:
    // конструкторы
    MATRIX()
    {
        n = m = 0;
        M = nullptr;
    }

    // конструктор с двумя параметрами
    MATRIX(int _m, int _n)
    {
        m = _m;
        n = _n;

        // Выделить память для матрицы
        // Выделить память для массива указателей
        M = (T**) new T * [m]; // количество строк, количество указателей

        // Выделить память для каждого указателя
        for (int i = 0; i < m; i++)
            M[i] = (T*)new T[n];

        // заполнить массив M нулями
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M[i][j] = 0;
    }

    // Конструктор копирования
    MATRIX(const MATRIX& _M)
    {
        // Создается новый объект для которого виделяется память
        // Копирование данных *this <= _M
        m = _M.m;
        n = _M.n;

        // Выделить память для M
        M = (T**) new T * [m]; // количество строк, количество указателей

        for (int i = 0; i < m; i++)
            M[i] = (T*) new T[n];

        // заполнить значениями
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M[i][j] = _M.M[i][j];
    }

    // методы доступа
    T GetMij(int i, int j)
    {
        if ((m > 0) && (n > 0))
            return M[i][j];
        else
            return 0;
    }
    int Getm()
    {
        return m;
    }

    int Getn()
    {
        return n;
    }

    void SetMij(int i, int j, T value)
    {
        if ((i < 0) || (i >= m))
            return;
        if ((j < 0) || (j >= n))
            return;
        M[i][j] = value;
    }

    // метод, выводящий матрицу
    void Print(const char* ObjName)
    {
        cout << "Object: " << ObjName << endl;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
                cout << fixed << scientific << setprecision(4) << setw(20) << setfill(' ') << M[i][j] << "\t";
            cout << endl;
        }
        cout << "---------------------" << endl << endl;
    }
    // метод печати
    void txt(const char* filename)
    {
        ofstream out;
        out.open(filename);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
                out << fixed << scientific << setprecision(4) << setw(20) << setfill(' ') << M[i][j] << "\t";
            out << endl;
        }
        out.close();
    }
//Конструктор, считывающий матрицу из файла txt
    MATRIX(int Row, int Col, string filename)
    {
        ifstream in(filename);

        if (in.is_open())
        {
            m = Row;
            n = Col;
            M = (T**) new T * [m];

            for (int i = 0; i < m; i++)
                M[i] = (T*) new T[n];


            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    in >> M[i][j];

            in.close();
        }
        else
        {
            cout << "File not found.";
            m = n = 0;
            M = nullptr;
        }

    }
//Удаляет все строки кроме первых "num" штук
    friend MATRIX ReduceRow(const MATRIX& _M, int num)
    {
        if (num < 0)
        {
            cout << "Negative Number!" << endl;
        }

        MATRIX Temp( num, _M.n);

        for (int i = 0; i < Temp.m; i++)
            for (int j = 0; j < Temp.n; j++)
                Temp.M[i][j] = _M.M[i][j];
        return Temp;
    }
    //Удаляет все столбцы кроме первых "num" штук
    friend MATRIX ReduceCol(const MATRIX& _M, int num)
    {
        if (num < 0)
        {
            cout << "Negative Number!" << endl;
        }

        MATRIX Temp(_M.m, num);

        for (int i = 0; i < Temp.m; i++)
            for (int j = 0; j < Temp.n; j++)
                Temp.M[i][j] = _M.M[i][j];

        return Temp;
    }
    // оператор присваивания 
    MATRIX operator=(const MATRIX& _M)
    {
        //Старая память освобождается, значения записываются в новую
        if (n > 0)
        {
            for (int i = 0; i < m; i++)
                delete[] M[i];
        }

        if (m > 0)
            delete[] M;

        m = _M.m;
        n = _M.n;

        M = (T**) new T * [m]; 
        for (int i = 0; i < m; i++)
            M[i] = (T*) new T[n];
        
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M[i][j] = _M.M[i][j];

        return *this;
    }
    // Деструктор - освобождает память, выделенную для матрицы
    ~MATRIX()
    {
        if (n > 0)
        {
            for (int i = 0; i < m; i++)
                delete[] M[i];
        }

        if (m > 0)
            delete[] M;
    }
    //Транспонирование
    friend MATRIX Transp(const MATRIX& _M)
    {
        MATRIX temp(_M.n, _M.m);
        for (int i = 0; i < _M.n; i++)
            for (int j = 0; j < _M.m; j++)
                temp.M[i][j] = _M.M[j][i];
        return temp;
    }
    //Перемножение матриц
    friend  MATRIX operator*(const MATRIX& _M, const MATRIX& _P)
    {
        if(_M.n!=_P.m)
           cout<<"Number of Columns of Matrix №1 is not equal to number of Rows of Matrix №2"<<endl;

        MATRIX temp(_M.m, _P.n);

        for (int i = 0; i < _M.m; i++)
            for (int j = 0; j < _P.n; j++)
                for (int k = 0; k < _M.n; k++)
                    temp.M[i][j] += _M.M[i][k] * _P.M[k][j];

        return(temp);
    }
    //Решение задачи на нахождение собственных чисел и векторов методом вращений Якоби
    friend MATRIX EigenSolver(const MATRIX& S, MATRIX* EigVal,double eps,int flag)
    {
        double maxVal = 1.001*eps;
        T phi;
        int i,j,r;
        int i0 = 1;
        int j0 = 1;
        int k = 0;
        int N = S.m;
        MATRIX Temp(S);//матрица собственных чисел(стремится к диагональной)
        MATRIX EigVec(N, N); //собственные вектора
        *EigVal = MATRIX(1, N);
        T a, b;

        //начало цикла while
        while (abs(maxVal) > eps)
        {
            maxVal = 0;
            for (i = 0; i < N; i++)
            {
                for (j = i + 1; j < N; j++)
                {
                    if (abs(Temp.M[i][j]) > abs(maxVal))
                    {
                        maxVal = Temp.M[i][j];
                        i0 = i;
                        j0 = j;
                    }
                }
            }

            phi = 0.5 * atan((2 * Temp.M[i0][j0]) / (Temp.M[i0][i0] - Temp.M[j0][j0]));
            
            if (k == 0)
            {
                for (i = 0; i < N; i++)
                    EigVec.M[i][i] = 1;

                EigVec.SetMij(i0, i0, cos(phi));
                EigVec.SetMij(i0, j0, -sin(phi));
                EigVec.SetMij(j0, i0, sin(phi));
                EigVec.SetMij(j0, j0, cos(phi));
            }
            else
            {
                for (i = 0; i < N; i++)
                {
                    a = EigVec.M[i][i0]; b = EigVec.M[i][j0];

                    EigVec.M[i][i0] = a * cos(phi) + b * sin(phi);
                    EigVec.M[i][j0] = -a * sin(phi) + b * cos(phi);
                }
            }
                                
            for (i = 0; i < N; i++)
            {
                a = Temp.M[i][i0]; b = Temp.M[i][j0];

                Temp.M[i][i0] = a * cos(phi) + b * sin(phi);
                Temp.M[i][j0] = -a * sin(phi) + b * cos(phi);
            }
            for (i = 0; i < N; i++)
            {
                a = Temp.M[i0][i]; b = Temp.M[j0][i];

                Temp.M[i0][i] = a * cos(phi) + b * sin(phi);
                Temp.M[j0][i] = -a * sin(phi) + b * cos(phi);
            }                                 
            k++;
        }

        if (k == 1)
            cout << "Matrix is diagonal or its non-diagonal elements' absolute value is less than " << eps <<endl
            <<"Try changing error value for better precision "<<endl;

        for ( i = 0; i < N; i++)
            EigVal->M[0][i] = Temp.M[i][i];

        //Пузырьковая cортировка по убыванию
        for ( i = 1; i < N; ++i)
        {
            for ( r = 0; r < N - i; r++)
            {
                if (EigVal->M[0][r] < EigVal->M[0][r + 1])
                {
                    // Обмен местами
                    T temp = EigVal->M[0][r];
                    EigVal->M[0][r] = EigVal->M[0][r + 1];
                    EigVal->M[0][r + 1] = temp;

                    for ( j = 0; j < N; j++)
                    {
                        T temp = EigVec.M[j][r];
                        EigVec.M[j][r] = EigVec.M[j][r + 1];
                        EigVec.M[j][r + 1] = temp;
                    }
                }
            }
        }
        if (flag == 1)
        {
            cout << "Number of iterations k = " << k << endl;
        }
        
        return EigVec;
    }

};// конец класса MATRIX

//Интерфейс Строителя
class Builder
{
public:
    virtual ~Builder() {};
    virtual void ProduceY( int, int, int, int,string) = 0;
    virtual void ProduceEigenVecAndValues(MATRIX<double>,double) = 0;
    virtual void ProduceSVD(double,int) = 0;
    virtual void ProduceNewData(vector <MATRIX<double> >* ,int) = 0;
    virtual vector <MATRIX<double> >* GetProduct() = 0;
};
//Конкретный строитель
class ConcreteBuilder : public Builder
{
private:
    vector <MATRIX<double> >* SVD; //0-Исходная матрица
public:                            //1-Правые вектора
    ConcreteBuilder()              //2-Собственные числа
    {                              //3-Левые вектора
        this->Reset();
    }
    ~ConcreteBuilder()
    {
        delete SVD;
    }

    void Reset()
    {
        this->SVD = new vector<MATRIX<double> >();
    }
    //Создание матрицы для дальнейшей работы
    void ProduceY( int Row, int Col, int k, int iftest, string filename)  
    {
        ofstream out;
        srand(static_cast<unsigned int>(time(0)));
        MATRIX<double> M(Row, Col);
        MATRIX<double> x1(1, Row);
        MATRIX<double> x2(1, Col);

        if (iftest == 1)
        {
            MATRIX<double> Y(Row, Col, filename);
            for (int i = 0; i < Row; i++)
                for (int j = 0; j < Col; j++)
                    M.SetMij(i, j, (Y.GetMij(i, j)));
        }
        else
        {
            int choice;
            int m, n;
            double noiseamp;
            choice = 0;
            cout << "Type 1 if you want 1.0 / sqrt(Row * Col)*cos(mx1+nx2) function" << endl;
            cout << "Type 2 if you want 1.0 / sqrt(Row * Col)*exp(-(x1+x2)^2) function" << endl;
            cout << "Type 3 if you want 1.0 / sqrt(Row * Col)*cos(mx1+nx2) function with phase noise" << endl;
            cout << "Type 4 if you want 1.0 / sqrt(Row * Col)*exp(-(x1+x2)^2) function with amplitude noise" << endl;

            while (choice != 1 && choice != 2 && choice != 3 && choice != 4)
            {
                cin >> choice;
                if (choice != 1 && choice != 2 && choice != 3 && choice != 4)
                    cout << "Such choice does not exist. Try again!" << endl;
            }
            switch(choice)
            {
            case 1:

                for (int i = 0; i < Row; i++)
                    x1.SetMij(0,i,-PI + 2.0 * i * PI / (Row - 1));
                for (int i = 0; i < Col; i++)
                    x2.SetMij(0, i, -PI + 2.0 * i * PI / (Col - 1));

                cout << "Type m and n" << endl;
                cin >> m >> n;
                for (int i = 0; i < Row; i++)
                    for (int j = 0; j < Col; j++)
                        M.SetMij(i, j, 1.0 / sqrt(Row * Col) *cos(m*x1.GetMij(0,i) + n*x2.GetMij(0,j)));
                break;
            case 2:
                
                for (int i = 0; i < Row; i++)
                    x1.SetMij(0, i, -1.0 + 2.0 * i * 1.0 / (Row - 1));
                for (int i = 0; i < Col; i++)
                    x2.SetMij(0, i, -1.0 + 2.0 * i * 1.0 / (Col - 1));

                for (int i = 0; i < Row; i++)
                    for (int j = 0; j < Col; j++)
                        M.SetMij(i, j, 1 / sqrt(Row * Col) * exp(-pow(x1.GetMij(0,i) + x2.GetMij(0,j), 2)));
                break;
            case 3:
                
                for (int i = 0; i < Row; i++)
                    x1.SetMij(0, i, -PI + 2.0 * i * PI / (Row - 1));
                for (int i = 0; i < Col; i++)
                    x2.SetMij(0, i, -PI + 2.0 * i * PI / (Col - 1));

                cout << "Type m , n and phase noise amplitude" << endl;
                cin >> m >> n>> noiseamp;
                for (int i = 0; i < Row; i++)
                    for (int j = 0; j < Col; j++)
                        M.SetMij(i, j, 1.0 / sqrt(Row * Col) * cos(m * x1.GetMij(0, i) + n * x2.GetMij(0, j) + noiseamp*getRandomNumber(-1,1)*PI));
                break;
            case 4:
                
                cout << "Type noise amplitude" << endl;
                cin >> noiseamp;
                for (int i = 0; i < Row; i++)
                    x1.SetMij(0, i, -1.0 + 2.0 * i * 1.0 / (Row - 1 ));
                for (int i = 0; i < Col; i++)
                    x2.SetMij(0, i, -1.0 + 2.0 * i * 1.0 / (Col - 1));

                for (int i = 0; i < Row; i++)
                    for (int j = 0; j < Col; j++)
                        M.SetMij(i, j, (1+noiseamp*getRandomNumber(-1, 1)) / sqrt(Row * Col) * exp(-pow(x1.GetMij(0, i) + x2.GetMij(0, j), 2)));
                break;
            }
            M.txt("Function.txt");
        }
        
        if (k == 1)
            this->SVD->push_back(M);
        else
            this->SVD->push_back(Transp(M));
    }
    //Нахождение собственных векторов
    void ProduceEigenVecAndValues(MATRIX<double> Y,double eps)
    {
        ofstream out;
        MATRIX<double>* Eig = new MATRIX<double>;
        MATRIX<double> M(Y.Getm(), Y.Getn());
        M = EigenSolver(Y, Eig,eps,1);
        M.txt ("Result.txt");
        Eig->txt("Eigenvalues.txt");
        this->SVD->push_back(M);
        this->SVD->push_back(*Eig);

        delete Eig;
    }

    void ProduceSVD(double eps,int flag)
    {
        ofstream out;
        int Col = SVD->at(0).Getn();
        int reduct = Col ;
        MATRIX<double>* Eig = new MATRIX<double>;
        MATRIX<double> M(Col, Col);

        M = EigenSolver(Transp(SVD->at(0)) * (SVD->at(0)), Eig,eps,flag);

        for (int i = 0; i < Eig->Getn(); i++)
        {
            if (abs(Eig->GetMij(0, i)) < 1E-15)
            {
                reduct = i;
		cout << endl<<endl;
		cout <<"--------------------WARNING----------------------------"<<endl;
                cout << "Machine zero eigenvalues were found"<<endl;
		cout <<"Such eigenvalues and their eigenvectors were removed" << endl;
		cout << "Resulting number of eigenvalues is "<< i <<endl<<endl;
                break;
            }
        }
        if (reduct != Col)
        {
            *Eig = ReduceCol(*Eig, reduct);
            M = ReduceCol(M, reduct);
        }

        MATRIX<double> invEig(Eig->Getn(), Eig->Getn());
        for (int i = 0; i < Eig->Getn(); i++)
        {
            Eig->SetMij(0, i, sqrt(Eig->GetMij(0, i)));
            invEig.SetMij(i, i, 1.0 / Eig->GetMij(0, i));
        }
            

        this->SVD->push_back(M);
        this->SVD->push_back(*Eig);
        this->SVD->push_back(SVD->at(0) * SVD->at(1) * invEig);//Y*U*invEig
       
        delete Eig;
    }

    void ProduceNewData(vector <MATRIX<double> >* SVD,int k)
    {
        int num;
        MATRIX<double>M(SVD->at(3).Getn(), SVD->at(1).Getm());
        MATRIX<double> S(SVD->at(2).Getn(), SVD->at(2).Getn());

        cout << "Resulting Eigenvalues are :" << endl;
        for(int i = 0; i < S.Getm(); i++)
        {
            S.SetMij(i, i, SVD->at(2).GetMij(0, i));
            cout << i + 1 << " ----->" << fixed << scientific << setprecision(4) << setw(15) << setfill(' ')<< S.GetMij(i, i) << endl;
        }
           
        cout << "Type number of members you want to keep" << endl;
        cin >> num;
        
        S = ReduceCol(S, num);
        S = ReduceRow(S, num);
        SVD->at(1) = ReduceCol(SVD->at(1), num);
        SVD->at(3) = ReduceCol(SVD->at(3), num);

        if (k == 1)
        {
            M = SVD->at(3) * S * Transp(SVD->at(1));
            SVD->at(1).txt("RightVectors.txt");
            SVD->at(3).txt("LeftVectors.txt");
        }
        else
        {
            M = Transp(SVD->at(3) * S * Transp(SVD->at(1)));
            SVD->at(1).txt("LeftVectors.txt");
            SVD->at(3).txt("RightVectors.txt");
        }
            
        SVD->at(2) = ReduceCol(SVD->at(2), num);
        SVD->at(2).txt("EigenValues.txt");
        M.txt("Result.txt");

    }

    vector <MATRIX<double> >* GetProduct()
    {
        vector <MATRIX<double> >* Result = this->SVD;
        this->Reset();
        return Result;
    }

};//end of ConcreteBuilder

class Director
{
private:
    Builder* builder;
public:
    void set_builder(Builder* builder)
    {
        this->builder = builder;
    }

    void Build_SVD()
    {
        int Row, Col,k,iftest;
        string filename;
        double eps;
        iftest = 123;
        cout << "Type 0 if you want to have a run with a test function " << endl;
        cout << "Type 1 if you want to read data from a file " << endl;

        while (iftest != 0 && iftest != 1 )
        {
            cin >> iftest;
            if (iftest != 0 && iftest != 1)
                cout << "Such choice does not exist. Try again!" << endl;
        }

        if (iftest == 1)
        {
            cout << "Type the name of the file that you want to read data from " << endl;
            cout << "Example: 123.txt" << endl;
            cin >> filename;
        }
        
        cout << "Type number of Rows " << endl;
        cin >> Row;
        cout << "Type number of Columns " << endl;
        cin >> Col;
        cout << "Set error value" << endl;
        cin >> eps;

        
        if (Col <= Row)
            k = 1;
        else
            k = 0;

        this->builder->ProduceY(Row, Col,k,iftest, filename);

        cout << "Waiting for results..." << endl << endl;

        this->builder->ProduceSVD(eps,1);
       
        vector <MATRIX<double> >* SVD = this->builder->GetProduct();

        this->builder->ProduceNewData(SVD,k);
        
    }

    void FindEigenVectorsAndValues()
    {
        string filename;
        int Row, Col;
        double eps;
        cout << "Type filename (example: 123.txt)" << endl;
        cin >> filename;
        cout << "Type number of Rows" << endl;
        cin >> Row;
        cout << "Type number of Columns" << endl;
        cin >> Col;
        cout << "Type error value" << endl;
        cin >> eps;

        MATRIX<double> Y(Row, Col, filename);

        this->builder->ProduceEigenVecAndValues(Y, eps);
       
        vector <MATRIX<double> >* SVD = this->builder->GetProduct();

        SVD->at(0).txt("Result.txt");
        SVD->at(1).txt("EigenValues.txt");
    }
};

void Client(Director& director)
{

    int choice=123;
  
    ConcreteBuilder* Builder = new ConcreteBuilder();
    director.set_builder(Builder);

    cout<<"Type 0 if you want to solve eigenvalue problem for a symmetric matrix"<<endl;
    cout<<"Type 1 if you want to perform Singular Value Decomposition (SVD)"<<endl;

    while(choice!=0 && choice!=1)
    {
      cin >> choice;
      if(choice!=0 && choice!=1)
	    cout<<"Such choice does not exist. Try again!"<<endl;
    }
    if(choice==0)
      director.FindEigenVectorsAndValues();
    if(choice==1)
      director.Build_SVD();

    delete Builder;
};

int main()
{
  //SimpleTimer timer;

    Director* director = new Director();

    Client(*director);

    delete director;

    return 0;
}
