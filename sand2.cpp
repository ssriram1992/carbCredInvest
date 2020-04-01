#include <iostream>

#include <map>

int main()
{
	std::map <int, int> Fibo;
	Fibo[1] = 1;
	Fibo[2] = 1;
	for(int i=3; i<10; i++) 
		Fibo[i] = Fibo[i-1] + Fibo[i-2];
	for (const auto x:Fibo){
		std::cout<<x.second<<" ";
	}
	std::cout<<'\n';
	return 0;
}
