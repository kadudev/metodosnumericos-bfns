#include <iostream> //biblioteca padrão do c++
#include <locale> //para uso da codificação em portugês
#include <cmath> //para uso da matemática mais "avançada no cód"
#include <iomanip> //para manipulação de casas decimais
#include <vector> //para uso de vetores para armazenamento de erros
#include <fstream> //para salvar erros em arquivo

 using namespace std;

 //função f(x) = x^2 -1
 double func(double x){
    return x*x -1;
 }

 double funcder(double x){
    return 2*x;
 }

 //bisseção
 double BissecaoMetodo(double xl, double xu, double eps, vector<double>& ErrosBissecao){
    
    if(func(xl)*func(xu) >=0){
        //se os limites do intervalo forem maior ou igual a zero, não há troca de sinais, e portanto a raiz não está nesse intervalo
        cout<<"Intervalo inválido. f(xl) * f(xu) deve ser menor que zero."<<endl;
        return NAN; 
    } 

    double xr; //novo xr, tal que xr = xl+xu/2
    double erro;
    int iteracao = 1;

    do{
        xr = (xl+xu)/2.0; 
        erro = fabs(func(xr)); //erro = |f(xr)|

        ErrosBissecao.push_back(erro); //salva o erro de cada iteração
        //imprimir a iteração e o erro atual
        cout<<"Iteração "<<iteracao<<": xr = "<<xr<<", f(xr) = "<<func(xr)<<", erro = "<<erro<<endl;

        //decide em qual lado do intervalo continua
        if(func(xr)*func(xl)<0){
            xu = xr; //o novo limite superior é o xr
        }
        else{
            xl = xr;
        }

        iteracao++; //incrementa o numero de iterações

    }while (erro > eps);
    return xr; //retorna a raiz aproximada
 }

//falsa posição
double FalsaPosicaoMetodo(double xl, double xu, double eps, vector<double>& ErrosFalsaPosicao){
    if(func(xl)*func(xu) >=0){
        //segue a mesma regra dos limites da bisseção
        cout<<"Intervalo inválido. xl * xu devem ser menor que zero."<<endl;
        return NAN; 
    } 

    double xr; //novo xr, tal que xr = xu - (f(xu)*(xl-xu)/(f(xl)-f(xu)))
    double erro;
    int iteracao = 1;

    do {
        //cálculo do xr
        xr = xu - (func(xu) * (xl - xu)) / (func(xl) - func(xu));
        erro = fabs(func(xr)); // |f(xr)|
        ErrosFalsaPosicao.push_back(erro); //salva o erro de cada iteração

        cout << "Iteração " << iteracao << ": xr = " << xr << ", f(xr) = " << func(xr) << ", erro = " << erro << endl;

        //atualiza os limites com base no sinal, da mesma forma que na bisseção
        if (func(xr) * func(xl) < 0) {
            xu = xr;
        } else {
            xl = xr;
        }

        iteracao++;

    } while (erro > eps);

    return xr;
}

//Método de Newton
double NewtonMetodo(double x0, double eps, vector<double>& ErrosNewton){
    double xi, erro;
    int iteracao = 1;
    do {
        if (funcder(x0) == 0) {
            cout << "Derivada é zero. Método de Newton falhou." << endl;
            return NAN;
        }

        xi = x0 - func(x0) / funcder(x0);
        erro = abs(xi - x0);
        x0 = xi;
        ErrosNewton.push_back(erro); //salva o erro de cada iteração
        cout << "Iteração " << iteracao << ": x"<<iteracao<<"= " << xi << ", f(x"<<iteracao<<") = " << func(xi) << ", erro = " << erro << endl;
        iteracao++;

    } while (erro > eps);

    return xi;
}

//método da secante
double SecanteMetodo(double x0, double x1, double eps, vector<double>& ErrosSecante){
    double x2, erro;
    int iteracao=1;
    cout<<"x0 = "<<x0<<", x1 = "<<x1<<endl;
    do{
        if(func(x1)-func(x0)==0){
            cout<<"Divisão por zero na secante"<<endl; //se f(x1) - f(x0) for zero, a divisão vai ser por zero
            return NAN;
        }

        x2=x1- func(x1)*(x1-x0) /(func(x1)-func(x0));
        erro = abs(x2-x1);
        x0=x1;
        x1=x2;
        ErrosSecante.push_back(erro); //salva o erro de cada iteração
        cout << "Iteração " << iteracao 
        << ": x"<<iteracao<<" = " << x2 
        << ", f(x"<<iteracao<<") = " << func(x2) 
        << ", erro = " << erro << endl;

        iteracao++;
    }while(erro>eps);
    return x2;   
}

void SalvarErrosCSV(const vector<double>& Bissecao, const vector<double>& FalsaPosicao, const vector<double>& Newton, const vector<double>& Secante) {
    //essa é uma função para salvar os erros em um arquivo de planilha, não faz parte da parte de cálculo dos métodos numéricos
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
    cout<<"Função utilizada para desenvolvimento dos métodos: f(x) = x² -1"<<endl;

    //declaraçãoo de variáveis para método da bisseção e falsa posição
    double xl = 0.0, xu = 2.0; //intervalo inicial [xl, xu]
    double eps = 1e-8; //precisão eps = 10^-8
    double x0=0.5; //chute inicial para o método de Newton e secante
    double x1=2.0; //segundo chute para o método da secante
    vector<double> ErrosBissecao; //vetor para guardar erros
    vector<double> ErrosFalsaPosicao; //vetor para guardar erros
    vector<double> ErrosNewton; //vetor para guardar erros
    vector<double> ErrosSecante; //vetor para guardar erros

    cout<<"---Método da bisseção---"<<endl<<"------------------------------------------------------"<<endl;
    double raizBissecao = BissecaoMetodo(xl, xu, eps, ErrosBissecao);
    cout << fixed << setprecision(12); //imprime valor com precisão de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (método da bisseção): "<<raizBissecao<<endl<<"------------------------------------------------------"<<endl<<endl;

    cout<<"---Método da falsa posição---"<<endl<<"------------------------------------------------------"<<endl;
    double raizFalsaPosicao = FalsaPosicaoMetodo(xl, xu, eps, ErrosFalsaPosicao);
    cout << fixed << setprecision(12); //imprime valor com precisão de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (método da falsa posição): "<<raizFalsaPosicao<<endl<<"------------------------------------------------------"<<endl<<endl;

    cout<<"---Método de Newton---"<<endl<<"------------------------------------------------------"<<endl;
    double raizNewton = NewtonMetodo(x0, eps, ErrosNewton);
    cout << fixed << setprecision(12); //imprime valor com precisão de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (método de Newton): "<<raizNewton<<endl<<"------------------------------------------------------"<<endl<<endl;

    cout<<"---Método da secante---"<<endl<<"------------------------------------------------------"<<endl;
    double raizSecante = SecanteMetodo(x0, x1, eps, ErrosSecante);
    cout << fixed << setprecision(12); //imprime valor com precisão de 12 casas decimais, para valor mais preciso
    cout<<"Raiz aproximada encontrada (método da secante): "<<raizSecante<<endl<<"------------------------------------------------------"<<endl<<endl;

    SalvarErrosCSV(ErrosBissecao, ErrosFalsaPosicao, ErrosNewton, ErrosSecante);
    return 0;
    
 }