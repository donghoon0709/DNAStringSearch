#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstring>
#include <ctime>

const char primerForward[3][50] = { "GACCCCAAAATCAGCGAAAT","TTACAAACATTGGCCGCAAA","AGATTTGGACCTGCGAGCG" };
const char primerReverse[3][50] = { "CAGATTCAACTGGCAGTAACCAGA","TTCTTCGGAATGTCGCGC","ACTTGTGGAGACAGCCGCTC" };

void NSS(char * sequence);
void KMP(char * sequence);

int main()
{
	FILE * dataFile = fopen("CoVID_Sequence.txt", "r");
	if (dataFile == NULL) {
		puts("File Not Found");
		return 0;
	}

	fseek(dataFile, 0, SEEK_END);
	int fileSize = ftell(dataFile);
	fseek(dataFile, 0, SEEK_SET);

	char * sequence = new char[fileSize];
	sequence[0] = '\0';
	char tmp[500];

	while (fgets(tmp, 100, dataFile) != NULL) {
		int tmpSize = strlen(tmp);
		tmp[tmpSize - 1] = '\0';
		strcat(sequence, tmp);
	}
	puts("-----Naive String Search 알고리즘-----");
	NSS(sequence);
	puts("--------------------------------------");

	puts("-----KMP 알고리즘-----");
	KMP(sequence);
	puts("----------------------");

	fclose(dataFile);
	delete sequence;

	return 0;
}

void NSS(char * sequence)
{
	clock_t start, end;
	int seqSize = strlen(sequence);

	start = clock();

	for (int primerIdx = 0; primerIdx < 3; ++primerIdx) {
		int primerForwardSize = strlen(primerForward[primerIdx]);
		int primerForwardCount = 0;
		int primerReverseSize = strlen(primerReverse[primerIdx]);
		int primerReverseCount = 0;

		for (int seq = 0; seq < seqSize - primerForwardSize; ++seq) {
			int primer;

			for (primer = 0; primer < primerForwardSize; ++primer) {
				if (sequence[seq + primer] != primerForward[primerIdx][primer]) break;
			}

			if (primer == primerForwardSize) {
				++primerForwardCount;
				printf("%d번째 염기에서 프라이머 %d의 Forward 발견, 프라이머 %d의 Forward 총 %d개\n", seq + 1, primerIdx + 1, primerIdx + 1, primerForwardCount);
			}
		}

		for (int seq = 0; seq < seqSize - primerReverseSize; ++seq) {
			int primer;

			for (primer = 0; primer < primerReverseSize; ++primer) {
				if (sequence[seq + primer] != primerReverse[primerIdx][primer]) break;
			}

			if (primer == primerReverseSize) {
				++primerReverseCount;
				printf("%d번째 염기에서 프라이머 %d의 Reverse 발견, 프라이머 %d의 Reverse 총 %d개\n", seq + 1, primerIdx + 1, primerIdx + 1, primerReverseCount);
			}
		}

		printf("결과 : 프라이머 %d의 Forward %d개, 프라이머 %d의 Reverse %d개\n\n", primerIdx + 1, primerForwardCount, primerIdx + 1, primerReverseCount);
	}

	end = clock();
	
	double elapsedTime = (double)(end - start) / CLOCKS_PER_SEC;

	printf("소요된 시간 : %lf초\n", elapsedTime);
}

void KMP(char * sequence)
{
	clock_t start, end;
	int seqSize = strlen(sequence);
	
	start = clock();

	for (int primerIdx = 0; primerIdx < 3; ++primerIdx) {
		int primerForwardSize = strlen(primerForward[primerIdx]);
		int primerForwardCount = 0;
		int primerForwardJump[50];

		int primerReverseSize = strlen(primerReverse[primerIdx]);
		int primerReverseCount = 0;
		int primerReverseJump[50];

		int head = -1, tail = 0;
		primerForwardJump[tail] = head;

		while (tail < primerForwardSize) {
			if (head == -1 || (head >= 0 && primerForward[primerIdx][head] == primerForward[primerIdx][tail])) {
				++head, ++tail;
				primerForwardJump[tail] = head;
			}
			else {
				head = primerForwardJump[head];
			}
		}
		int seq = 0, primer = 0;
		while (seq < seqSize) {
			if (primer == -1 || (primer >= 0 && sequence[seq] == primerForward[primerIdx][primer])) ++seq, ++primer;
			else if (sequence[seq] != primerForward[primerIdx][primer]) {
				primer = primerForwardJump[primer];
			}

			if (primer == primerForwardSize) {
				++primerForwardCount;
				printf("%d번째 염기에서 프라이머 %d의 Forward 발견, 프라이머 %d의 Forward 총 %d개\n", seq - primerForwardSize + 1, primerIdx + 1, primerIdx + 1, primerForwardCount);

				primer = primerForwardJump[primer];
			}
		}

		head = -1, tail = 0;
		primerReverseJump[tail] = head;

		while (tail < primerReverseSize) {
			if (head == -1 || (head >= 0 && primerReverse[primerIdx][head] == primerReverse[primerIdx][tail])) {
				++head, ++tail;
				primerReverseJump[tail] = head;
			}
			else {
				head = primerReverseJump[head];
			}
		}
		seq = 0, primer = 0;
		while (seq < seqSize) {
			if (primer == -1 || (primer >= 0 && sequence[seq] == primerReverse[primerIdx][primer])) ++seq, ++primer;
			else if (sequence[seq] != primerReverse[primerIdx][primer]) {
				primer = primerReverseJump[primer];
			}

			if (primer == primerReverseSize) {
				++primerReverseCount;
				printf("%d번째 염기에서 프라이머 %d의 Reverse 발견, 프라이머 %d의 Reverse 총 %d개\n", seq - primerReverseSize + 1, primerIdx + 1, primerIdx + 1, primerReverseCount);

				primer = primerReverseJump[primer];
			}
		}

		printf("결과 : 프라이머 %d의 Forward %d개, 프라이머 %d의 Reverse %d개\n\n", primerIdx + 1, primerForwardCount, primerIdx + 1, primerReverseCount);
	}

	end = clock();

	double elapsedTime = (double)(end - start) / CLOCKS_PER_SEC;

	printf("소요된 시간 : %lf초\n", elapsedTime);
}

// 프로그램 실행: <Ctrl+F5> 또는 [디버그] > [디버깅하지 않고 시작] 메뉴
// 프로그램 디버그: <F5> 키 또는 [디버그] > [디버깅 시작] 메뉴
