#include "Braun-Robbinson.hpp"




using namespace boost::numeric::ublas;

class Parametrs
{
    std::map<std::string, double> variables;

public:


double getValue(const std::string& paramName)
    {
        auto it = variables.find(paramName);
        if (it != variables.end())
        {
            return it->second;
        }
        
        throw std::invalid_argument("Неправильное имя переменной");
    }
    void setValue(const std::string& paramName, double value)
    {
        variables[paramName] = value;
    }

};


double core_func(Parametrs value)
{    
    double result = value.getValue("a")*(value.getValue("x")*value.getValue("x"))+value.getValue("b")*(value.getValue("y")*value.getValue("y"))+value.getValue("c")*value.getValue("x")*value.getValue("y")+value.getValue("d")*value.getValue("x")+value.getValue("e")*value.getValue("y");
    return result;    
};

void calculate_x(Parametrs *value){
    double result = -(value->getValue("c")*value->getValue("y")+value->getValue("d"))/(2*value->getValue("a"));
    value->setValue("x",result);
};
void calculate_y(Parametrs *value){
    double result = (value->getValue("c")*value->getValue("d")-2*value->getValue("a")*value->getValue("e"))/(4*value->getValue("a")*value->getValue("b")-(value->getValue("c")*value->getValue("c")));
    value->setValue("y",result);
};
int check_condition(Parametrs value){

    if(2*value.getValue("a")<0 && 2*value.getValue("b")>0){
        std::cout<<"Convex-concave!"<<std::endl;
        return true;
    }else{
        return false;
    }
}


double func_x_y(Parametrs value , std::string flag)
{
    double H_x = 2*value.getValue("a")*value.getValue("x")+value.getValue("c")*value.getValue("y")+value.getValue("d");
    double H_y = 2*value.getValue("b")*value.getValue("y")+value.getValue("c")+value.getValue("x")+value.getValue("e");

    if(flag=="x"){
        return H_x;
    }else if (flag=="y")
    {
        return H_y;
    }
    else
    {
        std::cout<<"Error, flag does not support"<<std::endl;
    }
    return -1;
}


matrix <double> Harr_init(Parametrs value, double size)
{
    matrix <double> H_arr(size+1,size+1);
    for (double i = 0; i <= size; ++i)
    {
        for (double j = 0; j <= size; ++j)
        {   
            value.setValue("x", i / size);
            value.setValue("y", j / size);
            H_arr(i,j)= round(1000.0 * core_func(value)) / 1000.0;
        }
    }

    return H_arr;
}




int main() {
    
    int count_eq = 0;
    int k=0;
    double N=2;
    Parametrs param;

    param.setValue("a",-0.5);
    param.setValue("b",(0.5/12.0));
    param.setValue("c",(10.0/3.0));
    param.setValue("d",(-2.0/3.0));
    param.setValue("e",(-4.0/3.0));

    

    calculate_y(&param);
    calculate_x(&param);
    


    std::cout<<param.getValue("x")<<std::endl;
    std::cout<<param.getValue("y")<<std::endl;
    
    
    
    int flag = check_condition(param);
        if(flag == true)
        {
           double sedlo =  core_func(param);
           std::cout<<"Sedlo"<<sedlo<<"\n"<<std::endl;
        }else{
        };
    matrix <double> eps_arr = Harr_init(param, N);
   
        for (int j = 0; j <= N; ++j)
        {   
           

            // std::cout<<boost::format{"%2%/%1%/%3%/%2%/%1%/%3%/%2%"}<<"{["<<eps_arr[0][j]<<"]["<<eps_arr[0+1][j]<<"]["<<eps_arr[0+2][j]<<" ]}"<<std::endl;
            std::cout<<boost::format{"%1$4.3f / %2$4.3f / %3$4.3f\n"}%eps_arr(0,j)%eps_arr(0+1,j)%eps_arr(0+2,j)<<std::endl;
        
        }
        std::cout<<"\n"<<std::endl;

    N=3;
    matrix<double> H_arr = Harr_init(param, N);

   
        for (int j = 0; j <= N; ++j)
        {   
           

            // std::cout<<boost::format{"%2%/%1%/%3%/%2%/%1%/%3%/%2%"}<<"{["<<eps_arr[0][j]<<"]["<<eps_arr[0+1][j]<<"]["<<eps_arr[0+2][j]<<" ]}"<<std::endl;
             std::cout<<boost::format{"%1$4.3f / %2$4.3f / %3$4.3f\n"}%H_arr(0,j)%H_arr(0+1,j)%H_arr(0+2,j)<<std::endl;
        }


braun_robinson(H_arr,3);

return 0;

}