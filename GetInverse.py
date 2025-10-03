# 행렬식을 이용해 역행렬 계산
def GetInverseByDeterminant(M):
    if (len(M) == 1) :  # 1x1은 직접 계산
        if M[0][0] != 0 : return [ [1 / M[0][0]] ]
        else :
            print("행렬식이 0입니다. 역행렬이 없습니다.") 
            return None  # 0 예외처리
    det = CalcDeterminant(M)    # 행렬식 계산

    # 행렬식이 0이면 역행렬 없음
    if det == 0:
        print("행렬식이 0입니다. 역행렬이 없습니다.")
        return None
    
    # 2x2는 직접 계산
    if (len(M) == 2):
        return [[M[1][1]/det, -1*M[0][1]/det],
                [-1*M[1][0]/det, M[0][0]/det]]

    cofactorMatrix = GetCofactorMatrix(M)   # 여인수행렬
    I = GetTranspose(cofactorMatrix)        # 여인수행렬의 전치

    # (1/det) * adjugate
    n = len(I)
    for i in range(n):
        for j in range(n):
            I[i][j] = I[i][j] / det     # 모든 성분에 행렬식 나눠주기
    return I    # 최종 역행렬 반환

def GetMinorOfMatrix(M, i, j):
    return [row[:j] + row[j+1:] for row in (M[:i] + M[i+1:])]

# 전치행렬 구하기
def GetTranspose(M):
    n, m = len(M), len(M)
    T = [[0] * n for _ in range(m)] # 내보낼 전치 행렬
    for i in range(n):
        for j in range(m):
            T[j][i] = M[i][j]
    return T

# 여인수행렬 구하기
def GetCofactorMatrix(M):
    n = len(M)
    cofactorMatrix = []
    for i in range(n):
        cofactorRow = []
        for j in range(n):
            cofactorRow.append(((-1)**(i+j)) * CalcDeterminant(GetMinorOfMatrix(M, i, j)))
        cofactorMatrix.append(cofactorRow)
    return cofactorMatrix

# 행렬식 계산
def CalcDeterminant(M):
    n = len(M)
    if n == 2:  # 2x2는 직접 계산
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    
    if n == 3:  # 3x3은 사루스 법칙 사용
        return (M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1]) - (M[0][2]*M[1][1]*M[2][0] + M[0][0]*M[1][2]*M[2][1] + M[0][1]*M[1][0]*M[2][2])
    
    det = 0
    for c in range(len(M)):
        det += ((-1)**c) * M[0][c] * CalcDeterminant(GetMinorOfMatrix(M, 0, c))     # 항상 0번 행을 기준으로 함
    return det

# 가우스-조던 소거법으로 역행렬 계산
def GetInverseByGaussJordan(M):
    M2 = [row[:] for row in M]  # 원본 유지를 위한 deep copy
    n = len(M2[0])

    # 역행렬이 될 단위행렬 I
    I = [[float(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        I[i][i] = 1
    
    for i in range(n):
        # 주 대각 성분에 접근했는데, 0이면 계산할 수 없음
        # 그러면 본인의 아래쪽 행들과 비교해서 계산할 수 있는 행이 있다면 해당 행과 교체
        if (M2[i][i] == 0):
            swap = -1
            for k in range(i+1, n):
                if (M2[k][i] != 0):
                    swap = k
                    break
            if (swap == -1):    # 바꿀 수 있는 행이 없음 : 역행렬이 없음
                print("Pivot이 0이며 대체행이 없습니다. 역행렬이 없습니다.")
                return None
            M2[i], M2[swap] = M2[swap], M2[i]
            I[i], I[swap] = I[swap], I[i]
        # 피벗을 기준으로 다른 행 계산
        for j in range(n):
            if (i == j): continue
            ratio = M2[j][i]/M2[i][i]   # 피벗(주 대각 성분 i, i)을 기준으로 본인 아랫쪽의 행과의 비율 계산
            for k in range(n): 
                M2[j][k] = M2[j][k] - ratio * M2[i][k]  # 계산한 비율을 토대로 모든 열 연산
                I[j][k] = I[j][k] - ratio * I[i][k]     # 역행렬이 될 단위 행렬에도 똑같은 비율을 적용하여 계산
    
    # 주 대각성분 1로 만들기
    for i in range(n):
        divisor = M2[i][i]
        for j in range(n):
            I[i][j] = I[i][j]/divisor   # 원래의 행렬은 굳이 계산할 필요가 없다.
    
    return I    # 역행렬 반환


# 행렬 출력
def PrintMatrix(M):
    n = len(M)
    if n == 1 : 
        # print("[ %7.3f ]" %M[0][0])   # 이렇게 하니까 부동소수점 오차로 인해 0이여야 할 것이 가끔 -0.000으로 출력됨
        print("[ %7.3f   ]" %CleanRound(M[0][0], 3))
        return

    print("┌ ", end = "")
    for i in range(n):
        print("%7.3f " %CleanRound(M[0][i], 3), end = "")
    print("  ┐")

    for x in range(1, n-1):
        print("│ ",  end = "")
        for y in range(n):
            print("%7.3f " %CleanRound(M[x][y], 3),  end = "")
        print("  │")

    print("└ ", end = "")
    for i in range(n):
        print("%7.3f " %CleanRound(M[n-1][i]),  end = "")
    print("  ┘")

# 부동 소수점 오차로 인해 단순 비교하면 같은 값으로 보여도 다르다고 나옴
def CompareMatrix(M, N):
    for i in range(len(M)):
        for j in range(len(N)):
            if CleanRound(M[i][j], 3) != CleanRound(N[i][j], 3):
                return False
    return True

# 부동 소수점 오차로 인해 가끔 -0이 출력되는 문제를 해결
def CleanRound(num, r=3):
    eps=1e-9
    num = round(num, r)
    if abs(num) < eps:
        return 0.0
    return num

# A * A^-1 = I 임을 활용해서 역함수가 참인지 확인한다.
def CheckInverseIsCorrect(M, Orig):
    n = len(M)
    Mult = MatrixMultiply(M, Orig)
    I = [[float(0) for _ in range(n)] for _ in range(n)]
    for i in range(len(I)):
        I[i][i] = 1
    if CompareMatrix(Mult, I) : print("역행렬이 올바릅니다.")
    else : print("역행렬이 올바르지 않습니다.")

# 행렬 곱셈
def MatrixMultiply(M, N):
    n = len(M)
    # 결과 행렬 초기화
    R = [[float(0) for _ in range(n)] for _ in range(n)]

    # 행렬 곱셈
    for i in range(n):
        for j in range(n):
            for k in range(n):
                R[i][j] += M[i][k] * N[k][j]
    return R

# x = A^-1 * b를 이용해 연립방정식 Ax = b 풀기
def GetSolutionOfLinearSystem(I, b):
    n = len(I)
    
    x =  [0.0 for _ in range(n)]
    for i in range(n):
        for j in range(n):
            x[i] += I[i][j] * b[j]
    return x

# 역행렬 계산 모드
def GetInverseMatrix():
    while True:
        while True: # 몇 차원 행렬로 할 것인지 입력
            try:
                n = int(input("역행렬을 계산할 정방 행렬의 차수 (-1 To Exit) : "))
                if (n < 1 and n != -1):
                    print("양의 정수로 입력해주세요.")
                    continue
                break
            except:
                print("양의 정수로 입력해주세요.")
        if (n == -1):
            break

        # n * n 행렬 생성 및 0으로 채워넣기
        M = [[float(0) for _ in range(n)] for _ in range(n)]

        inputNumbers = []
        while len(inputNumbers) < n * n:
            # 숫자 입력 (띄어쓰기 구분 가능)
            inputLine = input("숫자 입력 : ").split()
            numbersFromLine = []
            # 입력받은 줄 내에서 해석 가능 범위까지 해석
            try :
                for num in inputLine:
                    numbersFromLine.append(float(num))
            except:
                print("숫자를 제대로 입력해주세요.")

            # 해석한 숫자들을 하나씩 추가
            try:
                for num in numbersFromLine:
                    inputNumbers.append(num)
            except:
                print("숫자를 제대로 입력해주세요.")

        # 숫자 할당
        for i in range(n):
            for j in range(n):
                M[i][j] = inputNumbers[i * n + j]

        print("입력한 행렬 : ")
        PrintMatrix(M)

        # 행렬식을 이용하여 역행렬 계산
        Determ = GetInverseByDeterminant(M)
        if (Determ != None) :
            print("행렬식으로 계산한 역행렬 : ")
            PrintMatrix(Determ)
            CheckInverseIsCorrect(Determ, M)

        # 가우스-조던 소거법으로 역행렬 계산
        GJ = GetInverseByGaussJordan(M)
        if (GJ != None):
            print("가우스-조던 소거법으로 계산한 역행렬 : ")
            PrintMatrix(GJ)
            CheckInverseIsCorrect(Determ, M)

        if (Determ == None or GJ == None) : continue

        if (CompareMatrix(GJ, Determ)): print("두 방식의 결과가 같습니다.")
        else: print("! 두 방식의 결과가 다릅니다 !")

# 방정식 해 계산 모드
def SolveLinearSystem():
    while True:
        while True: # 몇 차원 행렬로 할 것인지 입력
            try:
                n = int(input("Ax = b 일 때 A의 차수 (-1 To Exit) : "))
                if (n < 1 and n != -1):
                    print("양의 정수로 입력해주세요.")
                    continue
                break
            except:
                print("양의 정수로 입력해주세요.")
        if (n == -1):
            break

        # n * n 행렬 생성 및 0으로 채워넣기
        A = [[float(0) for _ in range(n)] for _ in range(n)]

        inputNumbers = []
        while len(inputNumbers) < n * n:
            # 숫자 입력 (띄어쓰기 구분 가능)
            inputLine = input("A 숫자 입력 : ").split()
            numbersFromLine = []
            # 입력받은 줄 내에서 해석 가능 범위까지 해석
            try :
                for num in inputLine:
                    numbersFromLine.append(float(num))
            except:
                print("숫자를 제대로 입력해주세요.")

            # 해석한 숫자들을 하나씩 추가
            try:
                for num in numbersFromLine:
                    inputNumbers.append(num)
            except:
                print("숫자를 제대로 입력해주세요.")

        # 숫자 할당
        for i in range(n):
            for j in range(n):
                A[i][j] = inputNumbers[i * n + j]

        # n차 벡터 생성 및 0으로 채워넣기
        b = [float(0) for _ in range(n)]

        inputNumbers = []
        while len(inputNumbers) < n:
            # 숫자 입력 (띄어쓰기 구분 가능)
            inputLine = input("b 숫자 입력 : ").split()
            numbersFromLine = []
            # 입력받은 줄 내에서 해석 가능 범위까지 해석
            try :
                for num in inputLine:
                    numbersFromLine.append(float(num))
            except:
                print("숫자를 제대로 입력해주세요.")

            # 해석한 숫자들을 하나씩 추가
            try:
                for num in numbersFromLine:
                    inputNumbers.append(num)
            except:
                print("숫자를 제대로 입력해주세요.")

        # 숫자 할당 ( 굳이 이 방식을 사용하는 이유는, 입력이 더 많을 때 칸 수를 맞추기 위해서 )
        for i in range(n):
            b[i] = inputNumbers[i]

        I = GetInverseByDeterminant(A)          # 행렬식을 통한 역행렬 계산
        if I != None:
            x1 = GetSolutionOfLinearSystem(I, b)    # 방정식의 해 계산

            print("행렬식을 이용한 방정식의 해 : \n[ " , end="")
            for i in range(len(x1)):
                print("%7.3f " %CleanRound(x1[i], 3), end="")
            print("  ]")

        I = GetInverseByGaussJordan(A)          # 가우스-조던 소거법을 통한 역행렬 계산
        if I != None:
            x2 = GetSolutionOfLinearSystem(I, b)    # 방정식의 해 계산

            print("가우스-조던 소거법을 이용한 방정식의 해 : \n[ " , end="")
            for i in range(len(x2)):
                print("%7.3f " %CleanRound(x2[i], 3), end="")
            print("  ]")

            isSame = True
            for i in range(len(x1)):
                if CleanRound(x1[i], 3) != CleanRound(x2[i], 3):
                    print("역행렬이 서로 달라 해가 다릅니다!")
                    isSame = False
                    break
            if isSame: print("역행렬과 방정식의 해가 일치합니다!")


def main():
    while True:
        print("사용할 기능을 입력하세요.")
        print("1 : 역행렬 구하기")
        print("2 : 연립방정식 해 구하기")
        mode = int(input("모드 입력 (-1 To Exit) : "))

        getInverseMatrix = False
        solveLinearSystem = False

        if mode == 1 : getInverseMatrix = True
        elif mode == 2 : solveLinearSystem = True
        elif mode == -1 : break
        else : print("값을 제대로 입력해주세요!")

        if getInverseMatrix:
            GetInverseMatrix()
        if solveLinearSystem:
            SolveLinearSystem()

if __name__ == "__main__":
    main()