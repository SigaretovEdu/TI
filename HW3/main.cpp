#include <bits/chrono.h>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cmath>
#include <string>
#include <set>
#include <map>
#include <random>
#include <chrono>

void pr(std::vector<std::vector<int>>& matr);

void gen_Pmatrix(std::vector<std::vector<int>>& matr, size_t n, size_t k, int base) {
	matr.resize(k, std::vector<int>(n, 0));
	for(size_t i = 0; i < k; ++i) {
		matr[i][i] = 1;
	}
	std::srand(time(nullptr));
	for(size_t i = 0; i < matr.size(); ++i) {
		for(size_t j = k; j < matr[0].size(); ++j) {
			matr[i][j] = rand() % base;
		}
	}
}

void tr(std::vector<std::vector<int>>& matr) {
	std::vector<std::vector<int>> trMatr;
	trMatr.resize(matr[0].size(), std::vector<int>(matr.size()));
	for(size_t i = 0; i < matr.size(); ++i) {
		for(size_t j = 0; j < matr[0].size(); ++j) {
			trMatr[j][i] = matr[i][j];
		}
	}
	matr = trMatr;
}

void neg(std::vector<std::vector<int>>& matr, int base) {
	for(size_t i = 0; i < matr.size(); ++i)
		for(size_t j = 0; j < matr[0].size(); ++j)
			matr[i][j] = (base - matr[i][j]) % base;
}

void gen_Cmatr(std::vector<std::vector<int>>& matr, std::vector<std::vector<int>>& new_matr, size_t n, size_t k, int base) {
	new_matr.resize(k, std::vector<int>(n - k));
	for(size_t i = 0; i < new_matr.size(); ++i) {
		for(size_t j = 0; j < new_matr[0].size(); ++j) {
			new_matr[i][j] = matr[i][j + k];
		}
	}
	tr(new_matr);
	neg(new_matr, base);
	
	for(size_t i = 0; i < new_matr.size(); ++i)
		new_matr[i].resize(n, 0);
	for(size_t i = 0; i < n - k; ++i) {
		new_matr[i][i + k] = 1;
	}
}

void pr(std::vector<std::vector<int>>& matr) {
	for(size_t i = 0; i < matr.size(); ++i) {
		std::cout<<'\t';
		for(size_t j = 0; j < matr[0].size(); ++j)
			std::cout<<matr[i][j]<<' ';
		std::cout<<'\n';
	}
}

std::vector<std::vector<int>> mul(std::vector<std::vector<int>>& m1, std::vector<std::vector<int>>& m2, int base) {
	if(m1[0].size() != m2.size())
		throw std::runtime_error("~Can't calculate multiplication~");
	std::vector<std::vector<int>> ans;
	ans.resize(m1.size(), std::vector<int>(m2[0].size(), 0));
	for(size_t i = 0; i < ans.size(); ++i) {
		for(size_t j = 0; j < ans[0].size(); ++j) {
			for(size_t p = 0; p < m1[0].size(); ++p) {
				ans[i][j] = (ans[i][j] + m1[i][p] * m2[p][j]) % base;
			}
		}
	}
	return ans;
}

bool isEmpty(std::vector<std::vector<int>> matr) {
	for(size_t i = 0; i < matr.size(); ++i) {
		for(size_t j = 0; j < matr[0].size(); ++j) {
			if(matr[i][j] != 0) {
				return false;
			}
		}
	}
	return true;
}


void gen_combinations(int base, int size, std::vector<std::vector<int>>& comb) {
	std::vector<int> tmp;
	for(int i = 0; i < std::pow(base, size); ++i) {
		tmp.clear();
		int num = i;
		for(int j = 0; j < size; ++j) {
			tmp.insert(tmp.begin(), num % base);
			num /= base;
		}
		comb.push_back(tmp);
	}
}


template<typename T>
void info(const std::vector<T>& v) {
	for(size_t i = 0; i < v.size(); ++i)
		std::cout<<v[i];
	std::cout<<'\n';
}


template<typename T>
void info(const std::vector<std::vector<T>>& v) {
	for(size_t i = 0; i < v.size(); ++i) {
		for(size_t j = 0; j < v[0].size(); ++j) {
			std::cout<<v[i][j];
		}
		std::cout<<'\n';
	}
}


template<typename T>
void info(const std::vector<std::vector<std::vector<T>>>& v) {
	for(size_t i = 0; i < v.size(); ++i) {
		for(size_t j = 0; j < v[0].size(); ++j) {
			for(size_t k = 0; k < v[0][0].size(); ++k) {
				std::cout<<v[i][j][k];
			}
			std::cout<<' ';
		}
		std::cout<<'\n';
	}
}


template<typename T>
void info(const std::set<std::vector<T>>& st) {
	for(const auto& item: st) {
		for(size_t i = 0; i < item.size(); ++i)
			std::cout<<item[i];
		std::cout<<'\n';
	}
}


template<typename T>
void info(std::map<std::vector<T>, std::vector<T>>& mp) {
	for(const auto& item: mp) {
		for(size_t i = 0; i < item.first.size(); ++i)
			std::cout<<item.first[i];
		std::cout<<" = ";
		for(size_t i = 0; i < item.second.size(); ++i)
			std::cout<<item.second[i];
		std::cout<<'\n';
	}
}


template<typename T>
void to_file(const std::vector<T>& v, std::string name) {
	std::ofstream out(name);
	for(size_t i = 0; i < v.size(); ++i)
		out<<v[i];
	out<<'\n';
}


template<typename T>
void to_file(const std::vector<std::vector<T>>& v, std::string name) {
	std::ofstream out(name);
	for(size_t i = 0; i < v.size(); ++i) {
		for(size_t j = 0; j < v[0].size(); ++j) {
			out<<v[i][j];
		}
		out<<'\n';
	}
}


template<typename T>
void to_file(const std::vector<std::vector<std::vector<T>>>& v, std::string name) {
	std::ofstream out(name);
	for(size_t i = 0; i < v.size(); ++i) {
		for(size_t j = 0; j < v[0].size(); ++j) {
			for(size_t k = 0; k < v[0][0].size(); ++k) {
				out<<v[i][j][k];
			}
			out<<' ';
		}
		out<<'\n';
	}
}


template<typename T>
void to_file(const std::set<std::vector<T>>& st, std::string name) {
	std::ofstream out(name);
	for(const auto& item: st) {
		for(size_t i = 0; i < item.size(); ++i)
			out<<item[i];
		out<<'\n';
	}
}


template<typename T>
void to_file(std::map<std::vector<T>, std::vector<T>>& mp, std::string name) {
	std::ofstream out(name);
	for(const auto& item: mp) {
		for(size_t i = 0; i < item.first.size(); ++i)
			out<<item.first[i];
		out<<" = ";
		for(size_t i = 0; i < item.second.size(); ++i)
			out<<item.second[i];
		out<<'\n';
	}
}


std::vector<int> vec_sum(std::vector<int>& v1, std::vector<int>& v2, int base) {
	if(v1.size() != v2.size())
		throw std::runtime_error("~Unable to sum~");
	std::vector<int> ans(v1.size());
	for(size_t i = 0; i < ans.size(); ++i)
		ans[i] = (v1[i] + v2[i]) % base;
	return ans;
}

int weight(std::vector<int>& vec) {
	int ans = 0;
	for(size_t i = 0; i < vec.size(); ++i)
		if(vec[i] != 0)
			++ans;
	return ans;
}


void vecs_sum(std::vector<int>& place, std::vector<int>& comb, std::vector<std::vector<int>>& g, int base) {
	for(size_t i = 0; i < comb.size(); ++i) {
		for(int j = 0; j < comb[i]; ++j) {
			place = vec_sum(place, g[i], base);
		}
	}
}


int to_dec(std::vector<int>& vec, int base) {
	int ans = 0;
	for(size_t i = 0; i < vec.size(); ++i) {
		ans += vec[i] * std::pow(base, vec.size() - i);
	}
	return ans;
}


std::vector<int> to_vec(int num, int base, size_t size) {
	std::vector<int> ans;
	for(size_t i = 0; i < size; ++i) {
		ans.insert(ans.begin(), num % base);
		num /= base;
	}
	return ans;
}


void clear_row(std::vector<std::vector<int>>& row) {
	for(size_t i = 0; i < row.size(); ++i)
		for(size_t j = 0; j < row[0].size(); ++j)
			row[i][j] = 0;
}


void gen_lex_combinations(int base, int size, std::vector<std::vector<int>>& comb) {
	std::vector<std::vector<int>> tmp;
	gen_combinations(base, size, tmp);
	for(int i = 0; i <= size; ++i) {
		for(size_t j = 0; j < tmp.size(); ++j) {
			if(weight(tmp[j]) == i)
				comb.push_back(tmp[j]);
		}
	}
}


void gen_table(std::vector<std::vector<int>>& g, std::vector<std::vector<std::vector<int>>>& table, int base) {
	table.clear();
	size_t n = g[0].size(), k = g.size();
	std::vector<std::vector<int>> comb;
	gen_combinations(base, k, comb);
	std::vector<std::vector<int>> row(comb.size());
	for(size_t i = 0; i < row.size(); ++i) {
		row[i].resize(n, 0);
		vecs_sum(row[i], comb[i], g, base);
	}
	std::set<std::vector<int>> words;
	for(size_t i = 0; i < row.size(); ++i)
		words.insert(row[i]);
	table.push_back(row);
	clear_row(row);
	std::vector<std::vector<int>> alph;
	gen_lex_combinations(base, n, alph);
	for(auto& item: alph) {
		if(words.find(item) == words.end()) {
			row[0] = item;
			words.insert(item);
			for(size_t l = 1; l < row.size(); ++l) {
				row[l] = vec_sum(item, table[0][l], base);
				words.insert(row[l]);
			}
			table.push_back(row);
			clear_row(row);
		}
	}
}


std::vector<int> mul_vm(std::vector<int>& vec, std::vector<std::vector<int>>& matr, int base) {
	std::vector<int> ans(matr[0].size(), 0);
	for(size_t i = 0; i < matr[0].size(); ++i) {
		for(size_t j = 0; j < vec.size(); ++j) {
			ans[i] = (ans[i] + vec[j] * matr[j][i]) % base;
		}
	}
	return ans;
}


void gen_sindroms(std::vector<std::vector<std::vector<int>>>& table, std::map<std::vector<int>, std::vector<int>>& sindroms, std::vector<std::vector<int>>& hT, int base) {
	for(size_t i = 0; i < table.size(); ++i) {
		int min = 1000000;
		size_t pos = -1;
		bool f = true;
		for(size_t j = 0; j < table[i].size(); ++j) {
			int wei = weight(table[i][j]);
			if(wei < min) {
				pos = j;
				min = wei;
			}
			else if(wei == min) {
				f = false;
			}
		}
		if(f) {
			std::vector<int> sin = mul_vm(table[i][pos], hT, base);
			sindroms[sin] = table[i][pos];
		}
	}
}


void border() {
	std::cout<<"\n";
	for(int i = 0; i < 40; ++i)
		std::cout<<'~';
	std::cout<<"\n";
}


void channel(std::vector<int>& vec, int base) {
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr1(0, static_cast<int>(vec.size()) - 1);
    std::uniform_int_distribution<> distr2(0, base - 1);

    size_t pos = static_cast<size_t>(distr1(gen));
	int num = vec[pos];
	while(num == vec[pos])
		num = distr2(gen);
	vec[pos] = num;
}


std::vector<int> gen_word(size_t k, int base) {
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, base - 1);

	std::vector<int> word;
    for(size_t i = 0; i < k; ++i) {
        word.push_back(distr(gen));
	}

	return word;
}



void encode(std::vector<int>& word, std::vector<std::vector<int>>& g, int base) {
	word = mul_vm(word, g, base);
}


std::vector<int> find_error_dist(std::vector<int>& word, std::vector<std::vector<std::vector<int>>>& table) {
	std::vector<int> ans;
	for(size_t i = 0; i < table.size(); ++i) {
		for(size_t j = 0; j < table[i].size(); ++j) {
			if(table[i][j] == word) {
				ans = table[i][0];
			}
		}
	}
	return ans;
}


std::vector<int> find_error_sind(std::vector<int>& word, std::map<std::vector<int>, std::vector<int>>& sindroms, std::vector<std::vector<int>>& hT, int base) {
	std::vector<int> ans;
	std::vector<int> sin = mul_vm(word, hT, base);
	if(sindroms.find(sin) != sindroms.end()) {
		ans = sindroms[sin];
	}
	return ans;
}


std::vector<int> record_time(std::vector<int>& word, std::vector<std::vector<std::vector<int>>>& table, std::map<std::vector<int>, std::vector<int>>& sindroms, std::vector<std::vector<int>>& hT, int base) {

	auto start_dist = std::chrono::high_resolution_clock::now();
	std::vector<int> ans_dist = find_error_dist(word, table);
	auto stop_dist = std::chrono::high_resolution_clock::now();

	auto start_sind = std::chrono::high_resolution_clock::now();
	std::vector<int> ans_sind = find_error_sind(word, sindroms, hT, base);
	auto stop_sind = std::chrono::high_resolution_clock::now();

	//if(ans_sind != ans_dist) {
	//	std::cout<<"dist = ";
	//	info(ans_dist);
	//	std::cout<<"sind = ";
	//	info(ans_sind);
	//}
	//else if(!ans_sind.empty()) {
	//	std::cout<<"Error is = ";
	//	info(ans_sind);
	//	std::vector<int> recover = vec_sum(word, ans_sind, base);
	//	std::cout<<"Recovered vector is = ";
	//	info(recover);
	//}
	//else {
	//	std::cout<<"Unable to recover original word\n";
	//}
	//std::cout<<"execution time of dist method = "<<std::chrono::duration_cast<std::chrono::nanoseconds>(stop_dist - start_dist).count()<<"ns\n";
	//std::cout<<"execution time of sind method = "<<std::chrono::duration_cast<std::chrono::nanoseconds>(stop_sind - start_sind).count()<<"ns\n";
	
	std::vector<int> res(3, 0);
	res[0] = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_dist - start_dist).count();
	res[1] = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_sind - start_sind).count();
	if(ans_dist != ans_sind)
		res[2] = 1;
	return res;
}


void get_data(std::vector<std::vector<int>>& g, std::vector<std::vector<int>>& hT, std::vector<std::vector<std::vector<int>>>& table, std::map<std::vector<int>, std::vector<int>>& sindroms, int base) {
	std::vector<int> average(3, 0);
	int rounds = 100;
	for(int i = 0; i < rounds; ++i) {
		std::vector<int> word = gen_word(2, 2);
		encode(word, g, 2);
		channel(word, 2);
		std::vector<int> tmp = record_time(word, table, sindroms, hT, 2);
		average[0] += tmp[0];
		average[1] += tmp[1];
		average[2] += tmp[2];
	}
	std::cout<<"Для "<<rounds<<" экспериментов:\n";
	std::cout<<"Среднее время выполнения декодирования методом стандартных расположений = ";
	std::cout<<static_cast<float>(average[0]) / static_cast<float>(rounds)<<"нс\n";
	std::cout<<"Среднее время выполнения декодирования методом синдромов = ";
	std::cout<<static_cast<float>(average[1]) / static_cast<float>(rounds)<<"нс\n";
	std::cout<<"Количество раз, когда методом стандартных расположений найден ответ, а методом синдромов нет = ";
	std::cout<<average[2]<<'\n';
}


int main() {
	// Матрица G1 (4,2)-код над F2
	std::vector<std::vector<int>> g4_2, h4_2, hT4_2;
	std::vector<std::vector<std::vector<int>>> table4_2;
	std::map<std::vector<int>, std::vector<int>> sindroms4_2;
	gen_Pmatrix(g4_2, 4, 2, 2);
	gen_Cmatr(g4_2, h4_2, 4, 2, 2);
	hT4_2 = h4_2;
	tr(hT4_2);
	gen_table(g4_2, table4_2, 2);
	gen_sindroms(table4_2, sindroms4_2, hT4_2, 2);

	std::cout<<"G1 (4,2)-код над F2\n";
	std::cout<<"\nпорождающая матрица:\n\n";
	info(g4_2);
	std::cout<<"\nпроверочная матрица:\n\n";
	info(h4_2);
	std::cout<<"\nтранспонированная проверочная матрица:\n\n";
	info(hT4_2);
	std::cout<<"\nстандартное расположение:\n\n";
	info(table4_2);
	std::cout<<"\nтаблица синдромов (синдром = вектор ошибки):\n\n";
	info(sindroms4_2);
	std::cout<<'\n';
	get_data(g4_2, hT4_2, table4_2, sindroms4_2, 2);
	border();


	// Матрица G2 (15,11)-код над F2
	std::vector<std::vector<int>> g15_11, h15_11, hT15_11;
	std::vector<std::vector<std::vector<int>>> table15_11;
	std::map<std::vector<int>, std::vector<int>> sindroms15_11;
	gen_Pmatrix(g15_11, 15, 11, 2);
	gen_Cmatr(g15_11, h15_11, 15, 11, 2);
	hT15_11 = h15_11;
	tr(hT15_11);
	gen_table(g15_11, table15_11, 2);
	gen_sindroms(table15_11, sindroms15_11, hT15_11, 2);

	std::cout<<"G2 (15,11)-код над F2\n";
	std::cout<<"\nпорождающая матрица:\n\n";
	info(g15_11);
	std::cout<<"\nпроверочная матрица:\n\n";
	info(h15_11);
	std::cout<<"\nтранспонированная проверочная матрица:\n\n";
	info(hT15_11);
	std::cout<<"\nстандартное расположение в файле table15_11.txt\n\n";
	to_file(table15_11, "table15_11.txt");
	std::cout<<"\nтаблица синдромов (синдром = вектор ошибки):\n\n";
	info(sindroms15_11);
	std::cout<<'\n';
	get_data(g15_11, hT15_11, table15_11, sindroms15_11, 2);
	border();


	// Матрица G3 (6,4)-код над F3
	std::vector<std::vector<int>> g6_4, h6_4, hT6_4;
	std::vector<std::vector<std::vector<int>>> table6_4;
	std::map<std::vector<int>, std::vector<int>> sindroms6_4;
	gen_Pmatrix(g6_4, 6, 4, 3);
	gen_Cmatr(g6_4, h6_4, 6, 4, 3);
	hT6_4 = h6_4;
	tr(hT6_4);
	gen_table(g6_4, table6_4, 3);
	gen_sindroms(table6_4, sindroms6_4, hT6_4, 3);

	std::cout<<"G3 (6,4)-код над F3\n";
	std::cout<<"\nпорождающая матрица:\n\n";
	info(g6_4);
	std::cout<<"\nпроверочная матрица:\n\n";
	info(h6_4);
	std::cout<<"\nтранспонированная проверочная матрица:\n\n";
	info(hT6_4);
	std::cout<<"\nстандартное расположение в файле table6_4.txt\n\n";
	to_file(table6_4, "table6_4.txt");
	std::cout<<"\nтаблица синдромов (синдром = вектор ошибки):\n\n";
	info(sindroms6_4);
	std::cout<<'\n';
	get_data(g6_4, hT6_4, table6_4, sindroms6_4, 3);
	border();


	// Матрица G4 (8,4)-код над F3
	std::vector<std::vector<int>> g8_4, h8_4, hT8_4;
	std::vector<std::vector<std::vector<int>>> table8_4;
	std::map<std::vector<int>, std::vector<int>> sindroms8_4;
	gen_Pmatrix(g8_4, 8, 4, 3);
	gen_Cmatr(g8_4, h8_4, 8, 4, 3);
	hT8_4 = h8_4;
	tr(hT8_4);
	gen_table(g8_4, table8_4, 3);
	gen_sindroms(table8_4, sindroms8_4, hT8_4, 3);

	std::cout<<"G4 (8,4)-код над F3\n";
	std::cout<<"\nпорождающая матрица:\n\n";
	info(g8_4);
	std::cout<<"\nпроверочная матрица:\n\n";
	info(h8_4);
	std::cout<<"\nтранспонированная проверочная матрица:\n\n";
	info(hT8_4);
	std::cout<<"\nстандартное расположение в файле table8_4.txt\n\n";
	to_file(table8_4, "table8_4.txt");
	std::cout<<"\nтаблица синдромов (синдром = вектор ошибки):\n\n";
	info(sindroms8_4);
	std::cout<<'\n';
	get_data(g8_4, hT8_4, table8_4, sindroms8_4, 3);
}
