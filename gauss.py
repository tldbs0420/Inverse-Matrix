
def printMatrix(M):
    n = len(M[0])
    if n == 1 : 
        print("[ %7.3f ]" %M[0][0])
        return

    print("┌ ", end = "")
    for i in range(n):
        print("%7.3f " %M[0][i], end = "")
    print("  ┐")

    for x in range(1, n-1):
        print("│ ",  end = "")
        for y in range(n):
            print("%7.3f " %M[x][y],  end = "")
        print("  │")

    print("└ ", end = "")
    for i in range(n):
        print("%7.3f " %M[n-1][i],  end = "")
    print("  ┘")


def main():
    while True:
        try:
            n = int(input("역행렬을 계산할 정방 행렬의 차수 : "))
            break
        except:
            print("숫자를 제대로 입력해주세요.")

    M = [[float(0) for _ in range(n)] for _ in range(n)]

    inputNumbers = []
    while len(inputNumbers) < n * n:
        inputLine = input("숫자 입력 : ").split()
        numbersFromLine = []
        try :
            for t in inputLine:
                numbersFromLine.append(float(t))
        except:
            print("숫자를 제대로 입력해주세요.")
        try:
            for num in numbersFromLine:
                inputNumbers.append(num)
        except:
            print("숫자를 제대로 입력해주세요.")

    
    for i in range(n):
        for j in range(n):
            M[i][j] = inputNumbers[i * n + j]
    printMatrix(M)

if __name__ == "__main__":
    main()