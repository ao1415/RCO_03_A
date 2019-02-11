#include <iostream>

#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>

using namespace std;

class MilliSecTimer {
private:

	std::chrono::time_point<std::chrono::high_resolution_clock> s;
	unsigned long long int startCycle = 0;
	long long time = 0;

	const double CyclePerMilliSec = 2794000.0;

#ifndef _MSC_VER
	unsigned long long int getCycle() const {
		unsigned int low, high;
		__asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
		return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
	}
#endif // _MSC_VER

public:

	/// <summary>
	/// コンストラクタ
	/// </summary>
	MilliSecTimer() = default;
	/// <summary>
	/// コンストラクタ
	/// </summary>
	/// <param name="_time">設定時間(ミリ秒)</param>
	MilliSecTimer(const std::chrono::milliseconds& _time) noexcept { time = _time.count(); }

	/// <summary>
	/// 時間を設定する
	/// </summary>
	/// <param name="_time">設定時間(ミリ秒)</param>
	void set(const std::chrono::milliseconds& _time) noexcept { time = _time.count(); }

	/// <summary>
	/// タイマーを開始させる
	/// </summary>
	void start() noexcept {
#ifdef _MSC_VER
		s = std::chrono::high_resolution_clock::now();
#else
		startCycle = getCycle();
#endif // _MSC_VER
	}

	/// <summary>
	/// 設定時間経過したかを得る
	/// </summary>
	/// <returns>経過していれば true, それ以外は false</returns>
	inline const bool check() const noexcept {
#ifdef _MSC_VER
		const auto e = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count() >= time;
#else
		return (getCycle() - startCycle) / CyclePerMilliSec >= time;
#endif // _MSC_VER
	}

	/// <summary>
	/// 設定時間経過したかを得る
	/// </summary>
	/// <returns>経過していれば true, それ以外は false</returns>
	operator bool() const noexcept { return check(); }

	/// <summary>
	/// 経過時間を取得する(ミリ秒)
	/// </summary>
	/// <returns>計測時間(ミリ秒)</returns>
	inline const long long interval() const noexcept {
#ifdef _MSC_VER
		const auto e = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();
#else
		return static_cast<long long int>((getCycle() - startCycle) / CyclePerMilliSec);
#endif // _MSC_VER
	}

};

inline int readInt() {
	int val;
	cin >> val;
	return val;
}

struct Point {

	Point() {}
	Point(int x, int y) {
		this->x = x;
		this->y = y;
	}

	int x = -1;
	int y = -1;

	const string toString() const { return to_string(x) + " " + to_string(y); }

	Point operator+(const Point& o) const { return Point(x + o.x, y + o.y); }
	Point operator-(const Point& o) const { return Point(x - o.x, y - o.y); }

	void operator+=(const Point& o) { x += o.x; y += o.y; }
	void operator-=(const Point& o) { x -= o.x; y -= o.y; }

	const bool operator==(const Point& o) const { return (x == o.x && y == o.y); }
	const bool operator!=(const Point& o) const { return !(*this == o); }

	const bool operator<(const Point& o) const {
		if (y != o.y) return y < o.y;
		return x < o.x;
	}
};

struct XorShift {
	unsigned int x;
	XorShift() : x(2463534242U) {}
	unsigned int rand() {
		x ^= (x << 13);
		x ^= (x >> 17);
		x ^= (x << 5);
		return x;
	}

};

const inline double range(const Point& pos1, const Point& pos2) {
	const auto d = pos1 - pos2;

	return sqrt(d.x*d.x + d.y*d.y);
}

class AI {
private:

	const int limit = 1950;//ms

	const double TempStart = 10000.0;
	const double TempEnd = 0.01;
	const double Time = limit;
	const double TempDiff = (TempStart - TempEnd) / Time;

	const int N;
	const vector<Point>& points;

	XorShift random;

	bool probability(const double& base, const double& next, const long long& t) {

		const double diff = base - next;

		if (diff > 0) return true;

		const double temp = TempStart - TempDiff * t;

		const double p = exp(diff / temp) * 4294967295.0;

		return p > random.rand();
	}

	vector<double> getRanges(const vector<int>& answer) {

		const int size = (int)answer.size();
		vector<double> ranges(size, 0);

		for (int i = 0; i < size - 1; i++)
		{
			ranges[i] = range(points[answer[i]], points[answer[i + 1]]);
		}
		ranges[size - 1] = range(points[answer[size - 1]], points[answer[0]]);

		return ranges;
	}
	double getRanges(const vector<int>& answer, int p1, int p2) {

		const double r = range(points[answer[p1]], points[answer[p2]]);

		return r;
	}

	double getVariance(const vector<int>& answer, const vector<double>& ranges) {

		const int size = (int)answer.size();

		double sum = 0;
		for (int i = 0; i < size; i++)
		{
			sum += ranges[i];
		}

		const double average = sum / size;

		double sum2 = 0;
		for (int i = 0; i < size; i++)
		{
			const double sub = ranges[i] - average;
			sum2 += sub * sub;
		}

		return sum2 / size;
	}

	double getScore(const vector<int>& answer, const vector<double>& ranges) {

		const auto v = getVariance(answer, ranges);


		return (1000000.0 / (1.0 + v));
	}

	vector<int> init() {
		vector<int> answer;

		const double ave = 250.0 + (random.rand() % 10000) / 1000.0;

		answer.push_back(0);

		vector<int> flag(N, 0);
		flag[0] = 1;

		for (int i = 1; i < N; i++)
		{
			const Point pos = points[answer[i - 1]];

			int minJ = 0;
			double minR = 999999;
			for (int j = 1; j < N; j++)
			{
				if (flag[j] == 0)
				{
					double r = range(pos, points[j]) + (random.rand() % 10000) / 1000.0;

					if (minR > abs(r - ave))
					{
						minR = abs(r - ave);
						minJ = j;
					}
				}
			}

			flag[minJ] = 1;
			answer.push_back(minJ);

		}

		return answer;
	}

public:

	AI(const int _N, const vector<Point>& _points) : N(_N), points(_points) {

	}

	const vector<int> think() {


		MilliSecTimer timer;
		timer.set(chrono::milliseconds(limit));


		vector<int> bestAnswer = init();
		auto ranges = getRanges(bestAnswer);
		double bestScore = getScore(bestAnswer, ranges);

		long long int count = 0;

		timer.start();
		while (!timer.check())
		{
			count++;

			vector<int> answer = init();
			auto r = getRanges(answer);
			double score = getScore(answer, r);


			if (bestScore < score)
			{
				bestScore = score;
				bestAnswer = answer;
			}
		}

		cerr << "Score:" << bestScore << endl;
		cerr << "Count:" << count << endl;

		return bestAnswer;
	}

};

int main() {

	const int N = readInt();

	vector<Point> points(N);

	for (auto& pos : points)
	{
		pos.x = readInt();
		pos.y = readInt();
	}

	AI ai(N, points);

	const auto answer = ai.think();

	for (const auto& ans : answer)
	{
		cout << ans << endl;
	}

	return 0;
}
