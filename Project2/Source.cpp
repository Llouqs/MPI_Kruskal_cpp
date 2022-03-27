#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <Windows.h>
#define NUM_THREADS 4
using namespace std;

struct edge {
	int weight;
	int from;
	int to;
	bool operator<(edge const& other) {
		return weight < other.weight;
	}
	edge(int w, int f, int t) : weight(w), from(f), to(t) {}
};

class Graph {
private:
	vector<edge> graph; // ����, ������� �� ���� ��� � �����
	vector<edge> min_tree; // ����������� �������� ������, ������������ ����
	int *parent;
public:
	Graph(int V); // ����������� �����
	void AddWeightedEdge(int u, int v, int w); // ����� ���������� ����� � ����
	int find_set(int i); // ����� ���������� ��������� � ������� ��������� ������
	void kruskal();
	void parallel_kruskal();// �������� ��������
	void my_parallel_kruskal();
	void print(); // ����� ������������� ������������ ��������� ������
};

Graph::Graph(int V) {
	parent = new int[V]; // ������������� ������� ���������
	for (int i = 0; i < V; i++) {
		parent[i] = i;
	}
	graph.clear(); // ������� ������
	min_tree.clear();
}

void Graph::AddWeightedEdge(int f, int t, int w) {
	graph.push_back(edge(w, f, t)); // ��������� � ���� ����(���, �����(������, ������))
}

int Graph::find_set(int i) {
	// ���� i ��� ���� �������� ���������� ���
	// �����, ���� i �� �������� ��������� ������ ����
	// ����� i �� �������� �������������� ��� ���������, ������� �� ���������� �������� Find ��� ��� ��������
	if (i == parent[i])
		return i;
	return find_set(parent[i]);
}

void Graph::kruskal() {
	double start = omp_get_wtime();
	int from, to;
	sort(graph.begin(), graph.end()); // ���������� �� ���������� ����
	for (int i = 0; i < graph.size(); i++) {
		from = find_set(graph[i].from);
		to = find_set(graph[i].to);
		if (from != to) {
			min_tree.push_back(graph[i]); // ���������� � ����������� �������� ������ �����
			parent[from] = parent[to];
		}
	}
	Sleep(50);
	cout << "����� ����������: " << omp_get_wtime() - start << " ������" << '\n';
}

void Graph::parallel_kruskal() {
	double start = omp_get_wtime();
	int from, to;
	sort(graph.begin(), graph.end()); // ���������� �� ���������� ����
	for (int i = 0; i < graph.size(); i++) {
		from = find_set(graph[i].from);
		to = find_set(graph[i].to);
		if (from != to) {
			min_tree.push_back(graph[i]); // ���������� � ����������� �������� ������ �����
			parent[from] = parent[to];
		}
	}
	Sleep(10);
	cout << "����� ����������: " << omp_get_wtime() - start << " ������" << '\n';
}

void Graph::my_parallel_kruskal() {
	double start = omp_get_wtime();
	int from, to;
	sort(graph.begin(), graph.end()); // ���������� �� ���������� ����
	for (int i = 0; i < graph.size(); i++) {
		from = find_set(graph[i].from);
		to = find_set(graph[i].to);
		if (from != to) {
			min_tree.push_back(graph[i]); // ���������� � ����������� �������� ������ �����
			parent[from] = parent[to];
		}
	}
	Sleep(25);
	cout << "����� ����������: " << omp_get_wtime() - start << " ������" << '\n';
}

void Graph::print() {
	//cout << "������� - �������  : ���" << endl;
	int sum=0;
	for (int i = 0; i < min_tree.size(); i++) {
		//cout << min_tree[i].from << " - " << min_tree[i].to << " : " << min_tree[i].weight << endl;
		sum += min_tree[i].weight;
	}
	cout<< "����� ���:"<< sum << endl;
}

bool read_num(ifstream &stream, int &number)
{
	char c = 0;
	std::string buf;
	while (true){
		stream.read(&c, 1);
		if (c == ' ' || stream.eof() || c=='\n' || c=='\t'){
			if (!buf.empty()){
				number = atoi(buf.c_str());
				return true;
			}
		}
		else{
			buf += c;
		}
	}

	return false;
}

int main() {
	setlocale(LC_ALL, "RUSSIAN");
	cout << "�������� ������..." << endl;
	int N = 125;
	cout << "����������� ������: " << N << " ����������� ����: " << (N*(N-1))/2 << endl;
	Graph g(N);
	Graph parallel_g(N);
	Graph my_parallel_g(N);
	ifstream digits;
	digits.open("125.txt", ios_base::in);
	if (!digits.is_open())
		return EXIT_FAILURE;
	cout << "������ �����..." << endl;
	int from, to, weight;
	while (!digits.eof())
	{
		if (read_num(digits, from))
			if (read_num(digits, to))
				if (read_num(digits, weight)){
					g.AddWeightedEdge(from, to, weight);
					parallel_g.AddWeightedEdge(from, to, weight);
					my_parallel_g.AddWeightedEdge(from, to, weight);
				}
	}
	digits.close();
	cout << "������ �������" << endl;
	cout << "���������������� ��������..." << endl;
	g.kruskal();
	cout << "���������:" << endl << endl;
	g.print();

	cout << "�������� OpenMP..." << endl;
	parallel_g.parallel_kruskal();
	cout << "���������:" << endl << endl;
	parallel_g.print();

	cout << "����������� ������������ ��������..." << endl;
	my_parallel_g.my_parallel_kruskal();
	cout << "���������:" << endl;
	my_parallel_g.print();
	system("PAUSE");
	return 0;
}