#include <iostream>

using namespace std;

int main()
{
	cout<<"Enter a sequence of numbers to be summed:"<<endl;
	int num = 0,sum = 0;
	while(cin>>num)
	{
		if (num < 0) break;
		sum+=num;
	}
	cout<<"Sum is: "<<sum<<endl;
	
}