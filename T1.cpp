#include <iostream> //biblioteca padr�o do c++
#include <locale> //para uso da codifica��o em portug�s
#include <cmath> //para uso da matem�tica mais "avan�ada no c�d"
#include <iomanip> //para manipula��o de casas decimais
#include <vector> //para uso de vetores para armazenamento de erros
#include <fstream> //para salvar erros em arquivo

 using namespace std;

 //fun��o f(x) = x^2 -1
 double func(double x){
    return x*x -1;
 }

 double funcder(double x){
    return 2*x;
 }

 //bisse��o
 double BissecaoMetodo(double xl, double xu, double eps, vector<double>& ErrosBissecao){
    
    if(func(xl)*func(xu) >=0){
        //se os limites do intervalo forem maior ou igual a zero, n�o h� troca de sinais, e portanto a raiz n�o est� nesse intervalo
        cout<<"Intervalo inv�lido. f(xl) * f(xu) deve ser menor que zero."<<endl;
        return NAN; 
    } 

    double xr; //novo xr, tal que xr = xl+xu/2
    double erro;
    int iteracao = 1;

    do{
        xr = (xl+xu)/2.0; 
        erro = fabs(func(xr)); //erro = |f(xr)|

        ErrosBissecao.push_back(erro); //salva o erro de cada itera��o
        //imprimir a itera��o e o erro atual
        cout<<"Itera��o "<<iteracao<<": xr = "<<xr<<", f(xr) = "<<func(xr)<<", erro = "<<erro<<endl;

        //decide em qual lado do intervalo continua
        if(func(xr)*func(xl)<0){
            xu = xr; //o novo limite superior � o xr
        }
        else{
            xl = xr;
        }

        iteracao++; //incrementa o numero de itera��es

    }while (erro > eps);
    return xr; //retorna a raiz aproximada
 }

//falsa posi��o
double FalsaPosicaoMetodo(double xl, double xu, double eps, vector<double>& ErrosFalsaPosicao){
    if(func(xl)*func(xu) >=0){
        //segue a mesma regra dos limites da bisse��o
        cout<<"Intervalo inv�lido. xl * xu devem ser menor que zero."<<endl;
        return NAN; 
    } 

    double xr; //novo xr, tal que xr = xu - (f(xu)*(xl-xu)/(f(xl)-f(xu)))
    double erro;
    int iteracao = 1;

    do {
        //c�lculo do xr
        xr = xu - (func(xu) * (xl - xu)) / (func(xl) - func(xu));
        erro = fabs(func(xr)); // |f(xr)|
        ErrosFalsaPosicao.push_back(erro); //salva o erro de cada itera��o

        cout << "Itera��o " << iteracao << ": xr = " << xr << ", f(xr) = " << func(xr) << ", erro = " << erro << endl;

        //atualiza os limites com base no sinal, da mesma forma que na bisse��o
        if (func(xr) * func(xl) < 0) {
            xu = xr;
        } else {
            xl = xr;
        }

        iteracao++;

    } while (erro > eps);

    return xr;
}

//M�todo de Newton
double NewtonMetodo(double x0, double eps, vector<double>& ErrosNewton){
    double xi, erro;
    int iteracao = 1;
    do {
        if (funcder(x0) == 0) {
            cout << "Derivada � zero. M�todo de Newton falhou." << endl;
            return NAN;
        }

        xi = x0 - func(x0) / funcder(x0);
        erro = abs(xi - x0);
        x0 = xi;
        ErrosNewton.push_back(erro); //salva o erro de cada itera��o
        cout << "Itera��o " << iteracao << ": x"<<iteracao<<"= " << xi << ", f(x"<<iteracao<<") = " << func(xi) << ", erro = " << erro << endl;
        iteracao++;

    } while (erro > eps);

    return xi;
}

//m�todo da secante
double SecanteMetodo(double x0, double x1, double eps, vector<double>& ErrosSecante){
    double x2, erro;
    int iteracao=1;
    cout<<"x0 = "<<x0<<", x1 = "<<x1<<endl;
    do{
        if(func(x1)-func(x0)==0){
            cout<<"Divis�o por zero na secante"<<endl; //se f(x1) - f(x0) for zero, a divis�o vai ser por zero
            return NAN;
        }

        x2=x1- func(x1)*(x1-x0) /(func(x1)-func(x0));
        erro = abs(x2-x1);
        x0=x1;
        x1=x2;
        ErrosSecante.push_back(erro); //salva o erro de cada itera��o
        cout << "Itera��o " << iteracao 
        << ": x"<<iteracao<<" = " << x2 
        << ", f(x"<<iteracao<<") = " << func(x2) 
        << ", erro = " << erro << endl;

        iteracao++;
    }while(erro>eps);
    return x2;   
}

void SalvarErrosCSV(const vector<double>& Bissecao, const vector<double>& FalsaPosicao, const vector<double>& Newton, const vector<double>& Secante) {
    //essa � uma fun��o para salvar os erros em um arquivo de planilha, n�o faz parte da parte de c�lculo dos m�todos num�ricos
    ofstream arquivo("ErrosMetodos.csv");
    arquivo << "Iteracao,Bissecao,FalsaPosicao,Newton,Secante\n";

    size_t max_it = std::max(
        std::max(Bissecao.size(), FalsaPosicao.size()),
        std::max(Newton.size(), Secante.size())
    );
    for (size_t i = 0; i < max_it; ++i) {
        arquivo << i + 1 << ",";
        arquivo << (i < Bissecao.size() ? Bissecao[i] : 0) << ",";
        arquivo << (i < FalsaPosicao.size() ? FalsaPosicao[i] : 0) << ",";
        arquivo << (i < Newton.size() ? Newton[i] : 0) << ",";
        arquivo << (i < Secante.size() ? Secante[i] : 0) << "\n";
    }

    arquivo.close();
    cout << "Erros salvos em ErroMetodos.csv\n";
}


 int main()
 {
    setlocale(LC_ALL, "Portuguese");
    cout<<"Fun��o utilizada para desenvolvimento dos m�todos: f(x) = x� -1"<<endl;

    //declara��oo de vari�veis para m�todo da bisse��o e falsa posi��o
    double xl = 0.0, xu = 2.0; //intervalo inicial [xl, xu]
    double eps = 1e-8; //precis�o eps = 10^-8
    double x0=0.5; //chute inicial para o m�todo de Newton e secante
    double x1=2.0; //segundo chute para o m�todo da secante
    vector<double> ErrosBissecao; //vetor para guardar erros
    vector<double> ErrosFalsaPosicao; //vetor para guardar erros
    vector<double> ErrosNewton; //vetor para guardar erros
    vector<double> ErrosSecante; //vetor para guardar erros

    cout<<"---M�todo da bisse��o---"<<endl<<"------------------------------------------------------"<<endl;
    double raizBissecao = BissecaoMetodo(xl, xu, eps, ErrosBissecao);
    cout << fixed << setprecision(12); //imprime valor com precis�o de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (m�todo da bisse��o): "<<raizBissecao<<endl<<"------------------------------------------------------"<<endl<<endl;

    cout<<"---M�todo da falsa posi��o---"<<endl<<"------------------------------------------------------"<<endl;
    double raizFalsaPosicao = FalsaPosicaoMetodo(xl, xu, eps, ErrosFalsaPosicao);
    cout << fixed << setprecision(12); //imprime valor com precis�o de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (m�todo da falsa posi��o): "<<raizFalsaPosicao<<endl<<"------------------------------------------------------"<<endl<<endl;

    cout<<"---M�todo de Newton---"<<endl<<"------------------------------------------------------"<<endl;
    double raizNewton = NewtonMetodo(x0, eps, ErrosNewton);
    cout << fixed << setprecision(12); //imprime valor com precis�o de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (m�todo de Newton): "<<raizNewton<<endl<<"------------------------------------------------------"<<endl<<endl;

    cout<<"---M�todo da secante---"<<endl<<"------------------------------------------------------"<<endl;
    double raizSecante = SecanteMetodo(x0, x1, eps, ErrosSecante);
    cout << fixed << setprecision(12); //imprime valor com precis�o de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (m�todo da secante): "<<raizSecante<<endl<<"------------------------------------------------------"<<endl<<endl;

    SalvarErrosCSV(ErrosBissecao, ErrosFalsaPosicao, ErrosNewton, ErrosSecante);
    return 0;
    
 }