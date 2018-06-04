#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

enum Communication {
    PASS_UP_LEFT, PASS_UP, PASS_UP_RIGHT, PASS_RIGHT, PASS_DOWN_RIGHT, PASS_DOWN, PASS_DOWN_LEFT, PASS_LEFT
};

class Board {
public:
    Board(int _size, int _rank, int mpiSize);

    virtual ~Board();

    void clear();

    void nextRound();

    int readFile(const char *fileName);

    int writeFile(const char *fileName);

    int getRank();

    int getNodes();

    int getSize();

protected:
    int **board;
    int **nextBoard;
    int size;
    int nodes;
    int yStart;
    int yEnd;
    int xStart;
    int xEnd;
    int lowerLimit;
    int upperLimit;
    int rank;

    void setAlive(bool value, int x, int y);

    bool isAlive(int x, int y);

    int nextState(int i, int j);

    virtual void syncBorder() {};

    virtual void barrier() {};
};

class RowBoard : public Board {
public:
    RowBoard(int _size, int _rank, int mpiSize);

protected:
    virtual ~RowBoard() {};

    virtual void syncBorder();

};

class BlockBoard : public Board {
public:
    BlockBoard(int _size, int _rank, int mpiSize);

protected:
    int width;
    int height;
    int leftRank;
    int rightRank;

    virtual void syncBorder();

    void copyColumn(int *buffer, int index);

    void pasteColumn(int *buffer, int index);
};

BlockBoard::BlockBoard(int _size, int _rank, int mpiSize) : Board(_size, _rank, mpiSize) {
    int rowSquares = (int) sqrt(mpiSize);
    nodes = rowSquares * rowSquares;
    double subSquareSize = (double) _size / rowSquares;
    yStart = 1 + (rank % rowSquares) * subSquareSize;
    yEnd = 1 + ((rank % rowSquares) + 1) * subSquareSize;
    xStart = 1 + (rank / rowSquares) * subSquareSize;
    xEnd = 1 + ((rank / rowSquares) + 1) * subSquareSize;
    width = xEnd - xStart;
    height = yEnd - yStart;
    leftRank = _rank - rowSquares;
    rightRank = _rank + rowSquares;
    if (rank >= nodes) {
        xStart = xEnd = yStart = yEnd = lowerLimit = upperLimit = 0;
    }
}

Board::Board(int _size, int _rank, int mpiSize) {
    nodes = 1;
    size = _size + 2;
    yStart = lowerLimit = 1;
    yEnd = size - 1;
    upperLimit = yEnd - 1;
    board = new int *[size];
    nextBoard = new int *[size];
    rank = _rank;

    for (int i = 0; i < size; i++) {
        board[i] = new int[size];
        nextBoard[i] = new int[size];
        for (int j = 0; j < size; j++) {
            board[i][j] = 0;
            nextBoard[i][j] = 0;
        }
    }

    double subRangeLength = (double) _size / mpiSize;
    xStart = yStart;
    xEnd = yEnd;
    yStart = 1 + rank * subRangeLength;
    yEnd = 1 + (rank + 1) * subRangeLength;
}

Board::~Board() {
    for (int i = 0; i < size; i++) {
        delete[] board[i];
        delete[] nextBoard[i];
    }
    delete[] board;
    delete[] nextBoard;
}

int Board::getRank() {
    return rank;
}

int Board::getNodes() {
    return nodes;
}

int Board::getSize() {
    return size - 2;
}

void Board::setAlive(bool isAlive, int x, int y) {
    this->board[x][y] = isAlive ? 1 : 0;
}

bool Board::isAlive(int x, int y) {
    return board[x][y] == 1;
}

RowBoard::RowBoard(int _size, int _rank, int mpiSize) : Board(_size, _rank, mpiSize) {
    nodes = mpiSize;
};

void RowBoard::syncBorder() {
    MPI_Status statuses[4];
    int i = 0;
    if (rank % 2 == 0) {
        if (yStart > lowerLimit) {
            MPI_Send(board[yStart], size, MPI_INT, rank - 1, PASS_UP, MPI_COMM_WORLD);
            MPI_Recv(board[yStart - 1], size, MPI_INT, rank - 1, PASS_DOWN, MPI_COMM_WORLD, &statuses[i++]);
        }
        if (yEnd < upperLimit) {
            MPI_Send(board[yEnd - 1], size, MPI_INT, rank + 1, PASS_DOWN, MPI_COMM_WORLD);
            MPI_Recv(board[yEnd], size, MPI_INT, rank + 1, PASS_UP, MPI_COMM_WORLD, &statuses[i++]);
        }
    } else {
        if (yStart > lowerLimit) {
            MPI_Recv(board[yStart - 1], size, MPI_INT, rank - 1, PASS_DOWN, MPI_COMM_WORLD, &statuses[i++]);
            MPI_Send(board[yStart], size, MPI_INT, rank - 1, PASS_UP, MPI_COMM_WORLD);
        }
        if (yEnd < upperLimit) {
            MPI_Recv(board[yEnd], size, MPI_INT, rank + 1, PASS_UP, MPI_COMM_WORLD, &statuses[i++]);
            MPI_Send(board[yEnd - 1], size, MPI_INT, rank + 1, PASS_DOWN, MPI_COMM_WORLD);
        }
    }
}

void BlockBoard::copyColumn(int *buffer, int index) {
    for (int i = yStart, j = 0; i < yEnd; i++, j++) {
        buffer[j] = board[i][index];
    }
}

void BlockBoard::pasteColumn(int *buffer, int index) {
    for (int i = yStart, j = 0; i < yEnd; i++, j++) {
        board[i][index] = buffer[j];
    }
}

void BlockBoard::syncBorder() {
    MPI_Request requests[16];
    MPI_Status statuses[16];
    int i = 0;
    if (yStart > lowerLimit) {
        if (xStart > lowerLimit) {
            MPI_Isend(&board[yStart][xStart], 1, MPI_INT, leftRank - 1, PASS_UP_LEFT, MPI_COMM_WORLD, &requests[i++]);
            MPI_Irecv(&board[yStart - 1][xStart - 1], 1, MPI_INT, leftRank - 1, PASS_DOWN_RIGHT, MPI_COMM_WORLD,
                      &requests[i++]);
        }
        MPI_Isend(&board[yStart][xStart], width, MPI_INT, rank - 1, PASS_UP, MPI_COMM_WORLD, &requests[i++]);
        MPI_Irecv(&board[yStart - 1][xStart], width, MPI_INT, rank - 1, PASS_DOWN, MPI_COMM_WORLD, &requests[i++]);
        if (xEnd < upperLimit) {
            MPI_Isend(&board[yStart][xEnd - 1], 1, MPI_INT, rightRank - 1, PASS_UP_RIGHT, MPI_COMM_WORLD,
                      &requests[i++]);
            MPI_Irecv(&board[yStart - 1][xEnd], 1, MPI_INT, rightRank - 1, PASS_DOWN_LEFT, MPI_COMM_WORLD,
                      &requests[i++]);
        }
    }
    if (yEnd < upperLimit) {
        if (xStart > lowerLimit) {
            MPI_Isend(&board[yEnd - 1][xStart], 1, MPI_INT, leftRank + 1, PASS_DOWN_LEFT, MPI_COMM_WORLD,
                      &requests[i++]);
            MPI_Irecv(&board[yEnd][xStart - 1], 1, MPI_INT, leftRank + 1, PASS_UP_RIGHT, MPI_COMM_WORLD,
                      &requests[i++]);
        }
        MPI_Isend(&board[yEnd - 1][xStart], width, MPI_INT, rank + 1, PASS_DOWN, MPI_COMM_WORLD, &requests[i++]);
        MPI_Irecv(&board[yEnd][xStart], width, MPI_INT, rank + 1, PASS_UP, MPI_COMM_WORLD, &requests[i++]);
        if (xEnd < upperLimit) {
            MPI_Isend(&board[yEnd - 1][xEnd - 1], 1, MPI_INT, rightRank + 1, PASS_DOWN_RIGHT, MPI_COMM_WORLD,
                      &requests[i++]);
            MPI_Irecv(&board[yEnd][xEnd], 1, MPI_INT, rightRank + 1, PASS_UP_LEFT, MPI_COMM_WORLD, &requests[i++]);
        }
    }

    int outLeftBuffer[height];
    int inLeftBuffer[height];
    if (xStart > lowerLimit) {
        copyColumn(outLeftBuffer, xStart);
        MPI_Isend(outLeftBuffer, height, MPI_INT, leftRank, PASS_LEFT, MPI_COMM_WORLD, &requests[i++]);
        MPI_Irecv(inLeftBuffer, height, MPI_INT, leftRank, PASS_RIGHT, MPI_COMM_WORLD, &requests[i++]);
    }

    int outRightBuffer[height];
    int inRightBuffer[height];
    if (xEnd < upperLimit) {
        copyColumn(outRightBuffer, xEnd - 1);
        MPI_Isend(outRightBuffer, height, MPI_INT, rightRank, PASS_RIGHT, MPI_COMM_WORLD, &requests[i++]);
        MPI_Irecv(inRightBuffer, height, MPI_INT, rightRank, PASS_LEFT, MPI_COMM_WORLD, &requests[i++]);
    }

    MPI_Waitall(i, requests, statuses);
    if (xStart > lowerLimit) {
        pasteColumn(inLeftBuffer, xStart - 1);
    }
    if (xEnd < upperLimit) {
        pasteColumn(inRightBuffer, xEnd);
    }
}

int Board::nextState(int i, int j) {
    int sum = board[i - 1][j - 1] + board[i - 1][j] +
              board[i - 1][j + 1] + board[i][j - 1] +
              board[i][j + 1] + board[i + 1][j - 1] +
              board[i + 1][j] + board[i + 1][j + 1];

    if (sum < 2 || sum > 3) return 0;
    else if (sum == 3) return 1;
    else return board[i][j];
}

void Board::nextRound() {
    int i, j;
    syncBorder();

#pragma omp parallel private(i, j)
    {
#pragma omp for
        for (i = yStart; i < yEnd; i++) {
            for (j = xStart; j < xEnd; j++) {
                nextBoard[i][j] = nextState(i, j);
            }
        }
#pragma omp for
        for (int i = yStart; i < yEnd; i++) {
            for (int j = xStart; j < xEnd; j++) {
                board[i][j] = nextBoard[i][j];
            }
        }
    }
}

void Board::clear() {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            board[i][j] = 0;
            nextBoard[i][j] = 0;
        }
    }
}

int Board::readFile(const char *fileName) {
    ifstream file(fileName);
    char character;
    if (file.is_open()) {
        file >> noskipws;
        file >> character;
        for (int i = lowerLimit; i <= upperLimit && character != EOF; i++) {
            int j;
            for (j = lowerLimit; j <= upperLimit && character != '\n'; j++) {
                if (character == '1')
                    setAlive(true, i, j);
                else
                    setAlive(false, i, j);
                file >> character;
            }
            if (j >= size) {
                file.ignore(size, '\n');
            }
            file >> character;
        }
        file.close();
        return 0;
    } else {
        cout << "Unable to open file";
        return -1;
    }
}

int Board::writeFile(const char *fileName) {
    ofstream file(fileName);
    if (file.is_open()) {
        string prefix = "";
        for (int i = yStart; i < yEnd; i++) {
            file << prefix;
            for (int j = xStart; j < xEnd; j++) {
                if (isAlive(i, j))
                    file << "1";
                else
                    file << "0";
            }
            prefix = "\n";
        }
        file.close();
        return 0;
    } else {
        cout << "Unable to open file";
        return -1;
    }
}

void parseArgs(char **argv, int *size, int *iterations, int *threads) {
    if (argv[1] == NULL) {
        cout << "Insert board size [1st arg] " << endl;
        exit(-1);
    }
    sscanf(argv[1], "%d", size);

    if (argv[2] == NULL) {
        cout << "Insert iteration count [2nd arg] " << endl;
        exit(-1);
    }
    sscanf(argv[2], "%d", iterations);

    if (argv[3] == NULL) {
        cout << "Insert thread per node count [3rd arg] " << endl;
        exit(-1);
    }
    sscanf(argv[3], "%d", threads);

    if (argv[4] == NULL) {
        cout << "Insert input file name [4th arg] " << endl;
        exit(-1);
    }
}

void loadBoardFromFile(Board *board, char *infile) {
    board->clear();
    if (board->readFile(infile) < 0) {
        cout << "Missing or empty file input.txt" << endl;
        exit(-1);
    }
}

double runSingleThread(Board *board, int iterations, bool snapshotsOn) {
    double t1, t2, serialTime;
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    if (snapshotsOn) {
        int snapshots = 0;
        for (int i = 0; i < iterations; i++) {
            board->nextRound();
            snapshots++;
            std::stringstream outstream;
            outstream << "snapshot" << snapshots << ".txt";
            board->writeFile(outstream.str().c_str());
        }
    } else {
        t1 = omp_get_wtime();
        for (int i = 0; i < iterations; i++) {
            board->nextRound();
        }

        t2 = omp_get_wtime();
        serialTime = t2 - t1;
        printf("Single thread \nTime: %.3f\n", serialTime);
    }
    std::stringstream outstream;
    outstream << "output" << board->getSize() << ".txt";
    board->writeFile(outstream.str().c_str());
    return serialTime;
}

void runMultiThread(Board *board, int iterations, int threads, double serialTime, std::string label) {
    double t1, t2, time, speedUp, efficiency;
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = omp_get_wtime();
    for (int i = 0; i < iterations; i++) {
        board->nextRound();
    }
    t2 = omp_get_wtime();

    if (board->getRank() == 0) {
        time = t2 - t1;
        speedUp = serialTime / time;
        efficiency = speedUp / threads;
        printf("\n%s\n%d nodes * %d threads = %d threads \nTime: %.3f\nSpeedup: %f\nEfficiency: %f\n", label.c_str(),
               board->getNodes(), threads, threads * board->getNodes(), time, speedUp, efficiency);
    }
    std::stringstream outstream;
    outstream << "output" << board->getSize() << "-" << label.substr(1, 3) << board->getRank() << ".txt";
    board->writeFile(outstream.str().c_str());
}

int main(int argc, char *argv[]) {
    double serialTime;
    int size, iterations, threads = 1;
    int mpiSize, mpiRank, mpiTag, mpiStatus;
    bool snapshotsOn = argc == 6;
    char *infile = argv[4];

    mpiStatus = MPI_Init(&argc, &argv);
    if (mpiStatus != 0) {
        printf("Error while initializing MPI, status: %d", mpiStatus);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    parseArgs(argv, &size, &iterations, &threads);
    if (mpiRank == 0) {
        printf("The Game of Life\n");
        printf("N = %d\n", size);
        printf("Iterations = %d\n", iterations);
        printf("File name: %s\n", argv[4]);
        if (argv[5] == NULL) {
            cout << "Enter anything as a 5th argument to enable snapshotting" << endl << endl;
        }
        cout << "Starting the game..." << endl << endl;
        Board board = Board(size, mpiRank, 1);
        loadBoardFromFile(&board, infile);
        serialTime = runSingleThread(&board, iterations, false);
        if (snapshotsOn) {
            Board board = Board(size, mpiRank, 1);
            loadBoardFromFile(&board, infile);
            serialTime = runSingleThread(&board, iterations, snapshotsOn);
        }
    }

    Board *board = new RowBoard(size, mpiRank, mpiSize);
    loadBoardFromFile(board, infile);
    runMultiThread(board, iterations, threads, serialTime, "=ROW DECOMPOSITION=");
    delete board;
    board = new BlockBoard(size, mpiRank, mpiSize);
    loadBoardFromFile(board, infile);
    runMultiThread(board, iterations, threads, serialTime, "=BLK DECOMPOSITION=");
    delete board;
    MPI_Finalize();
    return 0;
}
