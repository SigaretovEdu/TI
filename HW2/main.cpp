#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>

void gen_vec(std::vector<int>& mes, size_t n, int q) {
	mes.clear();
	mes.resize(n, 0);
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, q - 1);

    for(size_t i = 0; i < mes.size(); ++i)
        mes[i] = distr(gen);
}

template<typename T>
void info(std::vector<T> v, int s = 1) {
	for(size_t i = 0; i < v.size(); ++i) {
		if(s == 1) {
			std::cout<<v[i]<<' ';
		}
		else {
			std::cout<<v[i]<<' ';
			for(int j = 0; j < s - 1; ++j) {
				std::cout<<"  ";
			}
		}
	}
	std::cout<<'\n';
}


void repeat_encoding(std::vector<int>& vec, size_t s) {
	std::vector<int> tmp;
	for(size_t i = 0; i < vec.size(); ++i) {
		for(size_t j = 0; j < s; ++j) {
			tmp.push_back(vec[i]);
		}
	}
	vec = tmp;
}


void repeat_decoding(std::vector<int>& vec, size_t s) {
	std::vector<int> tmp;
	for(size_t i = 0; i < vec.size(); i += s) {
		int pr = vec[i], cnt = 1;
		
		for(size_t pos = 1; pos < s; ++pos) {
			if(vec[i + pos] == pr) {
				++cnt;
			}
			else {
				--cnt;
				if(cnt == 0) {
					pr = vec[i + pos];
					cnt = 1;
				}
			}
		}

		tmp.push_back(pr);
	}
	vec = tmp;
}


void channel(std::vector<int>& vec, float p, int q) {
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> distr(0, q - 1);

	float coin;
	for(size_t i = 0; i < vec.size(); ++i) {
		coin = dis(gen);
		if(coin <= p) {
			int new_val = vec[i];
			while(new_val == vec[i]) {
        		new_val = distr(gen);
			}
			vec[i] = new_val;
		}
	}
}


int full_cycle(int n, int s, int q, float p) {
	std::vector<int> message;
	gen_vec(message, n, q);
	std::vector<int> copy = message;

	repeat_encoding(copy, s);
	
	channel(copy, p, q);
	
	repeat_decoding(copy, s);

	int errors = 0;
	for(size_t i = 0; i < message.size(); ++i)
		if(message[i] != copy[i])
			++errors;

	return errors;
}


float average(int n, int s, int q, float p) {
	float sum = 0.0;
	for(int i = 0; i < 10; ++i) {
		sum += full_cycle(n, s, q, p);
	}
	float avg = sum / 10.0;

	return avg;
}


void get_data1(std::vector<int>& x, std::vector<int>& y) {
	for(int i = 1; i <= 100; ++i) {
		x.push_back(i);
		y.push_back(average(i, 3, 5, 0.3));
	}
}


void get_data2(std::vector<float>& x, std::vector<int>& y) {
	for(float i = 0.0; i <= 1.0; i += 0.1) {
		x.push_back(i);
		y.push_back(average(100, 5, 5, i));
	}
}


void get_data3(std::vector<int>& x, std::vector<int>& y) {
	for(int i = 2; i <= 30; ++i) {
		x.push_back(i);
		y.push_back(average(100, i, 5, 0.3));
	}
}


int main() {
	/*
	 * Пункт 1, 2
	 * gen_vec(message, n, q) - генерирует сообщение длиной n, в q, записывает в вектор message
	 * repeat_encoding(message, s) - кодирование сообщения message кодом с повторениями с параметром s
	 * channel(message, p, q) - моделирует прохождение сообщения message через симметричный канал с параметром p
	 * repeat_decoding(message, s) - декодирование сообщения message кодом с повторениями с параметром s
	 * get_data1, get_data2, get_data3 - функции собирающие данные для построения графиков, для каждого эксперимернтального значения проводится 10 итераций и высичтывается среднее по количеству ошибок
	*/
	
	
	
	int n = 10;
	int s = 3;
	int q = 6;
	float p = 0.3;
	

	std::vector<int> message;
	gen_vec(message, n, q);
	std::cout<<"original message:\t";
	info(message, s);

	repeat_encoding(message, s);
	std::cout<<"encoded message:\t";
	info(message);
	
	channel(message, p, q);
	std::cout<<"noisy message:\t\t";
	info(message);
	
	repeat_decoding(message, s);
	std::cout<<"decoded message:\t";
	info(message, s);

	std::vector<int> x1, y1;
	get_data1(x1, y1);
	std::ofstream out1("../data1.txt");
	for(size_t i = 0; i < x1.size(); ++i)
		out1<<x1[i]<< ' ';
	out1<<'\n';
	for(size_t i = 0; i < y1.size(); ++i)
		out1 << y1[i] << ' ';
	out1.close();


	std::vector<float> x2;
	std::vector<int> y2;
	get_data2(x2, y2);
	std::ofstream out2("../data2.txt");
	for(size_t i = 0; i < x2.size(); ++i)
		out2<<x2[i]<< ' ';
	out2<<'\n';
	for(size_t i = 0; i < y2.size(); ++i)
		out2 << y2[i] << ' ';
	out2.close();

	std::vector<int> x3, y3;
	get_data3(x3, y3);
	std::ofstream out3("../data3.txt");
	for(size_t i = 0; i < x3.size(); ++i)
		out3<<x3[i]<< ' ';
	out3<<'\n';
	for(size_t i = 0; i < y3.size(); ++i)
		out3 << y3[i] << ' ';
	out3.close();
	
	return 0;
}
