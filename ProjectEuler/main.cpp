//
//  main.cpp
//  ProjectEuler
//
//  Created by ankit jain on 2/23/17.
//  Copyright © 2017 ankit jain. All rights reserved.
//

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <deque>
#include <stack>
#include "math.h"

using namespace std;

/*************************************************************************************************
*************************************************************************************************/

void Euler1()
{
    int sum = 0;
    for(int i=0; i<1000; i++)
    {
        if(i%3 == 0 || i%5 ==0)
            sum+=i;
    }
    
    cout<<sum;
}

/*************************************************************************************************
 *************************************************************************************************/

int fib(int n, vector<int> &greedy)
{
    if(greedy.size() > n)
    {
        return greedy[n];
    }
    
    int result = fib(n-1, greedy) + fib(n-2, greedy);
    greedy.push_back(result);
    
    cout<<result<<endl;
    return result;
}

void Euler2()
{
    // this will store all the previously calculated results
    vector<int> fibonacci;
    // add the first values
    fibonacci.push_back(0);
    fibonacci.push_back(1);
    
    int sum = 0;
    int k=2;
    while(1)
    {
        int result = fib(k, fibonacci);
        if(result > 4000000)
            break;
        else if(result%2==0)
            sum += result;
        
        k++;
    }
    
    cout<<sum<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

long primeFactors(long n, vector<long> &factors)
{
    // Print the number of 2s that divide n
    while (n%2 == 0)
    {
        factors.push_back(2);
        //printf("%d ", 2);
        n = n/2;
    }
    
    // n must be odd at this point.  So we can skip one element (Note i = i +2)
    for (long i = 3; i <= sqrt(n); i = i+2)
    {
        // While i divides n, print i and divide n
        while (n%i == 0)
        {
            factors.push_back(i);
            n = n/i;
        }
    }
    
    // This condition is to handle the case whien n is a prime number
    // greater than 2
    if(n>2)
        factors.push_back(n);
    
    return n;
}

void Euler3()
{
    vector<long> factors;
    cout<<primeFactors(600851475143, factors)<<endl;
    
    for(int i=0; i<factors.size(); i++)
        cout<<factors[i]<<endl;
}


/*************************************************************************************************
 *************************************************************************************************/

int largestPalindrome(int n)
{
    int upper_limit = 0;
    
    /* Loop to calculate upper bound(largest number
     of n-digit)*/
    for (int i=1; i<=n; i++)
    {
        upper_limit *= 10;
        upper_limit += 9;
    }
    
    // largest number of n-1 digit. One plus this number
    // is lower limit which is product of two numbers.
    int lower_limit = 1 + upper_limit/10;
    
    int max_product = 0; // Initialize result
    for (int i=upper_limit; i>=lower_limit; i--)
    {
        for (int j=i; j>=lower_limit; j--)
        {
            // calculating product of two n-digit numbers
            int product = i * j;
            if (product < max_product)
                break;
            int number = product;
            int reverse = 0;
            
            /* calculating reverse of product to check
             whether it is palindrome or not*/
            while (number != 0)
            {
                reverse = reverse * 10 + number % 10;
                number /= 10;
            }
            
            /* update new product if exist and if
             greater than previous one*/
            if (product == reverse && product > max_product)
                max_product = product;
        }
    }
    return max_product;
}

void Euler4()
{
    cout<<largestPalindrome(3)<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

int updateFactorMap(map<long, int> &factorMap, long factor, int count)
{
    // find if this factor already exists in the map
    map<long, int>::iterator it;
    it = factorMap.find(factor);
    
    // this is the updated difference in counts of this factor
    int diffCount = 0;
    
    if (it != factorMap.end())
    {
        if(it->second < count)
        {
            diffCount = count - it->second;
            it->second = count;
        }
    }
    else
    {
        factorMap.insert( pair<long, int>(factor, count) );
        diffCount = count;
    }
    
    return diffCount;
}

long updateFactorMap(map<long, int> &factorMap, long num)
{
    // find the prime factors of the number
    vector<long> factors;
    primeFactors(num, factors);
    
    // this is the updated multiple
    long updatedMultiple = 1;
    
    // once we have the prime factors update the map
    long factor = factors[0];
    int factorCount = 1;
    for(int i=1; i<factors.size(); i++)
    {
        if(factor == factors[i])
        {
            factorCount++;
        }
        else
        {
            int diffCount = updateFactorMap(factorMap, factor, factorCount);
            updatedMultiple *= powl(factor, diffCount);
            
            // update the factor and it's count
            factor = factors[i];
            factorCount = 1;
        }
    }
    
    int diffCount = updateFactorMap(factorMap, factor, factorCount);
    updatedMultiple *= powl(factor, diffCount);
    
    return updatedMultiple;
}

void Euler5()
{
    // max number
    long n=20;
    
    // smallest multiple
    long multiple = n;
    
    // a factor map which stores all the prime factors and count
    map<long, int> factorMap;
    updateFactorMap(factorMap, n);
    
    for(long i=n-1; i>1; i--)
    {
        // first check if the new number already divides the multiple
        if(multiple%i != 0)
        {
            multiple *= updateFactorMap(factorMap, i);
        }
    }
    
    cout<<multiple<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

long sumOfSquares(long n)
{
    return n*(n+1)*(2*n+1)/6;
}

long squareOfSum(long n)
{
    return n*n*(n+1)*(n+1)/4;
}

void Euler6()
{
    long val = 100;
    
    cout<<squareOfSum(val)-sumOfSquares(val)<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

bool testPrimality(long n)
{
    // Print the number of 2s that divide n
    while (n%2 == 0)
    {
        return false;
    }
    
    // n must be odd at this point.  So we can skip one element (Note i = i +2)
    for (long i = 3; i <= sqrt(n); i = i+2)
    {
        // While i divides n, print i and divide n
        if (n%i == 0)
        {
            return false;
        }
    }
    
    return true;
}

void Euler7()
{
    int count = 2;
    long i=6;
    long num = 0;
    
    while(1)
    {
        if( testPrimality(i-1) )
        {
            num = i-1;
            count++;
            
            if(count == 10001)
                break;
        }
        
        if( testPrimality(i+1) )
        {
            num = i+1;
            count++;
            
            if(count == 10001)
                break;
        }
            
        i+=6;
    }
    
    cout<<num<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler8()
{
    string s = "7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450";
    
    vector<int> numbers;
    for(int i=0; i<s.size(); i++)
    {
        numbers.push_back( s[i]-'0' );
    }
    
    int nLength = 13;
    long maxProduct = 1;
    
    // this queue will maintin a running product
    deque<long> runningProduct;
    
    // start iterating over the numbers
    for(int i=0; i<numbers.size(); i++)
    {
        // if a zero is encountered then purge the queue
        if(numbers[i] == 0)
        {
            runningProduct.clear();
            continue;
        }
        
        // if the running product has still not reached it's maximum size then do only push operation
        if(runningProduct.size() < nLength)
        {
            // start mainting running product in the queue
            if(runningProduct.empty())
                runningProduct.push_back(numbers[i]);
            else
                runningProduct.push_back(numbers[i]*runningProduct.back());
            
            // when queue reaches it's desired size
            if(runningProduct.size() == nLength)
            {
                if(maxProduct < runningProduct.back())
                    maxProduct = runningProduct.back();
            }
        }
        else
        {
            // get the latest product
            long currentProduct = runningProduct.back() * numbers[i] / numbers[i-nLength];
            
            if(currentProduct>maxProduct)
                maxProduct = currentProduct;
            
            // pop and push into the queue
            runningProduct.pop_front();
            runningProduct.push_back(currentProduct);
        }
    }
    
    cout<<maxProduct<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler9()
{
    int a=0, b=0, c=0;
    
    for(int i=1; i<500; i++)
    {
        if((1000*(500-i)) % (1000-i) == 0)
        {
            a = i;
            b = 1000*(500-a) / (1000-a);
            c = sqrtl(a*a+b*b);
        }
    }
    
    cout<<a*b*c<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler10()
{
    long sum = 2+3;
    
    for(long i=6; i<2000000; i+=6)
    {
        if(testPrimality(i-1))
            sum += i-1;
        
        if(testPrimality(i+1))
            sum += i+1;
    }
    
    cout<<sum<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler12()
{
    int numDivisors = 500;
    
    int triangleStart = 10;
    
    while(1)
    {
        long triangleNumber = triangleStart*(triangleStart+1)/2;
        
        // a factor map which stores all the prime factors and count
        map<long, int> factorMap;
        updateFactorMap(factorMap, triangleNumber);
        
        int count = 1;
        map<long, int>::iterator itr = factorMap.begin();
        
        for(; itr != factorMap.end(); ++itr)
        {
            count *= (itr->second+1);
        }
        
        if(count > numDivisors)
        {
            cout<<triangleNumber<<endl;
            return;
        }
        
        triangleStart++;
    }
}

/*************************************************************************************************
 *************************************************************************************************/

// this function is used to add an array of long numbers
vector<int> longAdd(vector<vector<int>> &add)
{
    int remainder = 0;
    int j=0;
    vector<int> digits;
    while(1)
    {
        // begin by initializing the sum of digits with previous remainder
        int sum = remainder;
        bool bFound = false;
        for(int i=0; i<add.size(); i++)
        {
            // calculate sum of all digits in one location
            if(j < add[i].size())
            {
                sum += add[i][j];
                bFound = true;
            }
        }
        
        if(bFound)
        {
            digits.push_back(sum%10);
            remainder = sum/10;
            
            j++;
        }
        else
            break;
    }
    
    while(remainder)
    {
        digits.push_back(remainder%10);
        remainder /= 10;
    }
    
    return digits;
}

void Euler13()
{
    string s[] =
    {
    "37107287533902102798797998220837590246510135740250",
    "46376937677490009712648124896970078050417018260538",
    "74324986199524741059474233309513058123726617309629",
    "91942213363574161572522430563301811072406154908250",
    "23067588207539346171171980310421047513778063246676",
    "89261670696623633820136378418383684178734361726757",
    "28112879812849979408065481931592621691275889832738",
    "44274228917432520321923589422876796487670272189318",
    "47451445736001306439091167216856844588711603153276",
    "70386486105843025439939619828917593665686757934951",
    "62176457141856560629502157223196586755079324193331",
    "64906352462741904929101432445813822663347944758178",
    "92575867718337217661963751590579239728245598838407",
    "58203565325359399008402633568948830189458628227828",
    "80181199384826282014278194139940567587151170094390",
    "35398664372827112653829987240784473053190104293586",
    "86515506006295864861532075273371959191420517255829",
    "71693888707715466499115593487603532921714970056938",
    "54370070576826684624621495650076471787294438377604",
    "53282654108756828443191190634694037855217779295145",
    "36123272525000296071075082563815656710885258350721",
    "45876576172410976447339110607218265236877223636045",
    "17423706905851860660448207621209813287860733969412",
    "81142660418086830619328460811191061556940512689692",
    "51934325451728388641918047049293215058642563049483",
    "62467221648435076201727918039944693004732956340691",
    "15732444386908125794514089057706229429197107928209",
    "55037687525678773091862540744969844508330393682126",
    "18336384825330154686196124348767681297534375946515",
    "80386287592878490201521685554828717201219257766954",
    "78182833757993103614740356856449095527097864797581",
    "16726320100436897842553539920931837441497806860984",
    "48403098129077791799088218795327364475675590848030",
    "87086987551392711854517078544161852424320693150332",
    "59959406895756536782107074926966537676326235447210",
    "69793950679652694742597709739166693763042633987085",
    "41052684708299085211399427365734116182760315001271",
    "65378607361501080857009149939512557028198746004375",
    "35829035317434717326932123578154982629742552737307",
    "94953759765105305946966067683156574377167401875275",
    "88902802571733229619176668713819931811048770190271",
    "25267680276078003013678680992525463401061632866526",
    "36270218540497705585629946580636237993140746255962",
    "24074486908231174977792365466257246923322810917141",
    "91430288197103288597806669760892938638285025333403",
    "34413065578016127815921815005561868836468420090470",
    "23053081172816430487623791969842487255036638784583",
    "11487696932154902810424020138335124462181441773470",
    "63783299490636259666498587618221225225512486764533",
    "67720186971698544312419572409913959008952310058822",
    "95548255300263520781532296796249481641953868218774",
    "76085327132285723110424803456124867697064507995236",
    "37774242535411291684276865538926205024910326572967",
    "23701913275725675285653248258265463092207058596522",
    "29798860272258331913126375147341994889534765745501",
    "18495701454879288984856827726077713721403798879715",
    "38298203783031473527721580348144513491373226651381",
    "34829543829199918180278916522431027392251122869539",
    "40957953066405232632538044100059654939159879593635",
    "29746152185502371307642255121183693803580388584903",
    "41698116222072977186158236678424689157993532961922",
    "62467957194401269043877107275048102390895523597457",
    "23189706772547915061505504953922979530901129967519",
    "86188088225875314529584099251203829009407770775672",
    "11306739708304724483816533873502340845647058077308",
    "82959174767140363198008187129011875491310547126581",
    "97623331044818386269515456334926366572897563400500",
    "42846280183517070527831839425882145521227251250327",
    "55121603546981200581762165212827652751691296897789",
    "32238195734329339946437501907836945765883352399886",
    "75506164965184775180738168837861091527357929701337",
    "62177842752192623401942399639168044983993173312731",
    "32924185707147349566916674687634660915035914677504",
    "99518671430235219628894890102423325116913619626622",
    "73267460800591547471830798392868535206946944540724",
    "76841822524674417161514036427982273348055556214818",
    "97142617910342598647204516893989422179826088076852",
    "87783646182799346313767754307809363333018982642090",
    "10848802521674670883215120185883543223812876952786",
    "71329612474782464538636993009049310363619763878039",
    "62184073572399794223406235393808339651327408011116",
    "66627891981488087797941876876144230030984490851411",
    "60661826293682836764744779239180335110989069790714",
    "85786944089552990653640447425576083659976645795096",
    "66024396409905389607120198219976047599490197230297",
    "64913982680032973156037120041377903785566085089252",
    "16730939319872750275468906903707539413042652315011",
    "94809377245048795150954100921645863754710598436791",
    "78639167021187492431995700641917969777599028300699",
    "15368713711936614952811305876380278410754449733078",
    "40789923115535562561142322423255033685442488917353",
    "44889911501440648020369068063960672322193204149535",
    "41503128880339536053299340368006977710650566631954",
    "81234880673210146739058568557934581403627822703280",
    "82616570773948327592232845941706525094512325230608",
    "22918802058777319719839450180888072429661980811197",
    "77158542502016545090413245809786882778948721859617",
    "72107838435069186155435662884062257473692284509516",
    "20849603980134001723930671666823555245252804609722",
    "53503534226472524250874054075591789781264330331690",
    };
    
    int nStrings = 100;
    int nLength = 50;
    
    // prepare an array of vectors each containing the long number
    // each long number itself is an array of it's digits
    vector<vector<int>> vecNumbers;
    
    for(int j=0; j<nStrings; j++)
    {
        vector<int> vecDigits;
        for(int i=nLength-1; i>-1; i--)
        {
            vecDigits.push_back(s[j][i] - '0');
        }
        
        vecNumbers.push_back(vecDigits);
    }
    
    vector<int> sumDigits = longAdd(vecNumbers);
    
    for(int i=(int)sumDigits.size()-1; i>sumDigits.size()-11; i--)
    {
        cout<<sumDigits[i];
    }
    cout<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

long paths(int row, int col, int &maxRow, int &maxCol, long** dp)
{
    if(dp[row][col])
        return dp[row][col];
    
    long val = 0;
    
    if(row==maxRow)
    {
        val = paths(row, col+1, maxRow, maxCol, dp);
        dp[row][col] = val;
        return val;
    }
    
    if(col==maxCol)
    {
        val = paths(row+1, col, maxRow, maxCol, dp);
        dp[row][col] = val;
        return val;
    }

    val = paths(row+1, col, maxRow, maxCol, dp) + paths(row, col+1, maxRow, maxCol, dp);
    dp[row][col] = val;
    return val;
}

void Euler15()
{
    int max = 20;
    
    long** dpPaths = new long*[max+1];
    for(int i=0; i<max+1; i++)
    {
        dpPaths[i] = new long[max+1];
    }
    
    for(int i=0; i<max+1; i++)
    {
        for(int j=0; j<max+1; j++)
        {
            dpPaths[i][j] = 0;
        }
    }
    
    dpPaths[max][max] = 0;
    dpPaths[max][max-1] = 1;
    dpPaths[max-1][max] = 1;
    
    cout<<paths(0, 0, max, max, dpPaths)<<endl;
    
    for(int i=0; i<max+1; i++)
    {
        delete [] dpPaths[i];
    }
    delete [] dpPaths;
}

/*************************************************************************************************
 *************************************************************************************************/

// create an array of digits from a number
vector<int> longNumber(int num)
{
    vector<int> digits;
    
    while(num)
    {
        digits.push_back(num%10);
        num /= 10;
    }
    
    return digits;
}

// this function is used to multiply one long number with a single digit
// if values of nZeroes is set then that many zeroes are inserted first in the number
vector<int> longMult(vector<int> &a, int b, int nZeroes = 0)
{
    vector<int> result;
    for (int i=0; i<nZeroes; i++) {
        result.push_back(0);
    }
    
    int carry = 0;
    for(int i=0; i<a.size(); i++)
    {
        int mult = a[i]*b+carry;
        result.push_back(mult%10);
        carry = mult/10;
    }
    
    while(carry)
    {
        result.push_back(carry%10);
        carry /= 10;
    }
    
    return result;
}

// this function is used to multiply two long numbers
vector<int> longMult(vector<int> &a, vector<int> &b)
{
    vector<vector<int>> vecNumbers;
    for(int i=0; i<b.size(); i++)
    {
        vecNumbers.push_back(longMult(a, b[i], i));
    }
    
    return longAdd(vecNumbers);
}

// this function takes in a long number and raises it to a power
vector<int> longPow(vector<int> &number, int p)
{
    vector<int> result(number);
    
    for(int i=0; i<p-1; i++)
    {
        result = longMult(result, number);
    }
    
    return result;
}

void Euler16()
{
    vector<int> digits = longNumber(256);
    vector<int> result(digits);
    result = longPow(result, 5);
    result = longPow(result, 5);
    result = longPow(result, 5);
    
    int sum=0;
    for(int i=0; i<result.size(); i++)
    {
        sum += result[i];
    }
    
    cout<<sum<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler17()
{
    int nDigits[] = {
        3,  // one
        3,  // two
        5,  // three
        4,  // four
        4,  // five
        3,  // six
        5,  // seven
        5,  // eight
        4   // nine
    };
    
    int nElevenToTwenty[] = {
        6,  // eleven
        6,  // twelve
        8,  // thirteen
        8,  // fourteen
        7,  // fifteen
        7,  // sixteen
        9,  // seventeen
        8,  // eighteen
        8   // nineteen
    };
    
    int nTens[] = {
        3,  // ten
        6,  // twenty
        6,  // thirty
        5,  // forty
        5,  // fifty
        5,  // sixty
        7,  // seventy
        6,  // eighty
        6   // ninety
    };
    
    int nHundred = 7;   // hundred
    int nThousand = 8;  // thousand
    int nAnd = 3;       // and
    
    int digitsLength = 0;
    
    // add first 9 digits (multiplied by ten) for each occurance in a decade
    for(int i=0; i<9; i++)
        digitsLength += 9*nDigits[i];
    
    // add Eleven to twenty
    for(int i=0; i<9; i++)
        digitsLength += nElevenToTwenty[i];
    
    // add 2nd digit (multiplied by nine for number of occurances
    for(int i=1; i<9; i++)
        digitsLength += 9*nTens[i];
    
    // add the tens
    for(int i=0; i<9; i++)
        digitsLength += nTens[i];
    
    int digitsLessThanHundred = digitsLength;
    
    // add the terms for hundreds
    for(int i=0; i<9; i++)
        digitsLength += 100*nDigits[i] + 100*nHundred + digitsLessThanHundred;
    
    // adding the word and
    digitsLength += (900-9)*nAnd;
    
    digitsLength += nThousand+nDigits[0];
    cout<<digitsLength<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

// returns maximum in an array
int maxArray(int *arr, int len)
{
    int max = arr[0];
    
    for(int i=1; i<len; i++)
        if(arr[i] > max)
            max = arr[i];
    
    return max;
}

int EvaluateTriangleSum(int **input, int **output, int row, int col)
{
    if(col < 0 || col > row || row < 0)
        return 0;
    
    if(output[row][col] == 0)
    {
        output[row][col] =  input[row][col] + std::max( EvaluateTriangleSum(input, output, row-1, col), EvaluateTriangleSum(input, output, row-1, col-1) );
    }
    
    return output[row][col];
}

void Euler18()
{
    FILE* fp = NULL;
    fp = fopen("triangle18.txt", "r");
    
    int nLines = 15;
    int **nTriangles = new int*[nLines];
    for(int i=0; i<nLines; i++)
    {
        nTriangles[i] = new int[i+1];
        for(int j=0; j<i+1; j++)
            fscanf(fp, "%d", &nTriangles[i][j]);
    }
    
    // this will be a recursive solution with DP
    // allocate memory for holding results
    int **nTriangleSum = new int*[nLines];
    for(int i=0; i<nLines; i++)
    {
        nTriangleSum[i] = new int[i+1];
        for(int j=0; j<i+1; j++)
            nTriangleSum[i][j] = 0;
    }
    
    for(int i=0; i<nLines; i++)
        EvaluateTriangleSum(nTriangles, nTriangleSum, nLines-1, i);
    
    cout<<maxArray(nTriangleSum[nLines-1], nLines)<<endl;
    
    for(int i=0; i<nLines; i++)
    {
        delete [] nTriangles[i];
        delete [] nTriangleSum[i];
    }
    
    delete [] nTriangles;
    delete [] nTriangleSum;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler20()
{
    vector<int> digits;
    digits.push_back(1);
    for(int i=2; i<101; i++)
    {
        vector<int> num = longNumber(i);
        digits = longMult(digits, num);
    }
    
    int sum = 0;
    for(int i=0; i<digits.size(); i++)
        sum += digits[i];
    
    cout<<sum<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

int sumOfDivisors(int num)
{
    // initialize sum
    int sum = 1;
    int sqrtNum = sqrt(num);
    
    // number till which loop will run
    int max = sqrtNum+1;
    
    // add the square root if perfect square
    if(sqrtNum*sqrtNum == num)
    {
        sum += sqrtNum;
        max--;
    }
    
    for(int i=2; i<max; i++)
    {
        if(num % i == 0)
            sum += i + num/i;
    }
    
    return sum;
}

void Euler21()
{
    // this is the number uptil which amicable numbers have to be summed
    int maxNum = 10000;
    
    // create a simple array map of numbers and their sum of divisors
    int *mapSumDivisiors = new int[maxNum];
    memset(mapSumDivisiors, 0, sizeof(int)*maxNum);
    
    // initialize desired sum
    int sumAmicable = 0;
    
    // fill in starting values in sum divisors
    mapSumDivisiors[0] = 0;
    mapSumDivisiors[1] = 0;
    mapSumDivisiors[2] = 1;
    mapSumDivisiors[3] = 1;
    
    // now iterate for the rest of the numbers
    for(int i=4; i<maxNum; i++)
    {
        int sumDivs = sumOfDivisors(i);
        mapSumDivisiors[i] = sumDivs;
        
        if(sumDivs < i && mapSumDivisiors[sumDivs] == i)
            sumAmicable += i+sumDivs;
    }
    
    cout<<sumAmicable<<endl;
    
    delete [] mapSumDivisiors;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler22()
{
    FILE* fp = NULL;
    fp = fopen("p022_names.txt", "r");
    
    char cc;
    // vector of all names in file
    vector<string> names;
    
    // boolean to track begin or end of name
    string name = "";
    while(cc!=EOF)
    {
        cc = getc(fp);
        if(cc == '"')
        {
            if(!name.empty())
            {
                names.push_back(name);
                name = "";
            }
        }
        else if(cc != ',')
        {
            name += cc;
        }
    }
    
    // sort the names
    sort(names.begin(), names.end());
    
    long sum = 0;
    
    for(int i=0; i<names.size(); i++)
    {
        int nameSum = 0;
        for(int j=0; j<names[i].length(); j++)
        {
            nameSum += 1 + names[i][j] - 'A';
        }
        
        sum += nameSum*(i+1);
    }
    
    cout<<sum<<endl;
    
    fclose(fp);
}

/*************************************************************************************************
 *************************************************************************************************/

// input is the number of which we want to get the divisor
// second input is a map which contains corrpesponding to each number it's set of divisors
// one big assumption is that we have already found out all the divisors of all the numbers less than the requested number
const unordered_set<int>& getDivisors(int num, unordered_map<int, unordered_set<int>> &divisorsMap)
{
    // check if this number already exists in the map
    unordered_map<int, unordered_set<int>>::const_iterator result = divisorsMap.find(num);
    
    // if found then return result as it is
    if(result != divisorsMap.end())
        return result->second;
    
    // get the sqrt of number
    int sqrtNum = sqrt(num);
    
    // number till which loop will run
    int max = sqrtNum+1;
    
    // check if this is a perfect square
    bool bPerfectSquare = false;
    if(sqrtNum*sqrtNum == num)
    {
        bPerfectSquare = true;
        max--;
    }
    
    for(int i=2; i<max; i++)
    {
        if(num % i == 0)
        {
            // retrieve the divisors of second number i.e num/i
            unordered_map<int, unordered_set<int>>::const_iterator divMap = divisorsMap.find(num/i);
            
            if(divMap == divisorsMap.end())
            {
                // not possible
                throw "divisors map not created correctly, need to debug";
            }
            
            // create a set from the divisors of second number
            unordered_set<int> divisors(divMap->second);
            
            // now walk over the previous divisors and update if necessary
            unordered_set<int>::const_iterator divItr = divMap->second.begin();
            for(; divItr != divMap->second.end(); ++divItr)
            {
                divisors.insert( (*divItr) * i );
            }
            
            // for good measure add the perfect square root, might not have been added if it was a prime
            if(bPerfectSquare)
                divisors.insert(sqrtNum);
            
            // for good measure add num/i as well, might not have been added if it was a prime
            divisors.insert(num/i);
            
            // add this in the map of divisors
            divisorsMap.insert(make_pair(num, divisors));

            // break now, since we must have updated the divisors
            return divisorsMap[num];
        }
    }
    
    // if reached this point then it means no divisor found yet
    unordered_set<int> divisors;
    divisors.insert(1);
    
    // that means number is prime or square of a prime
    if(bPerfectSquare)
        divisors.insert(sqrtNum);
    
    // add this in the map of divisors
    divisorsMap.insert(make_pair(num, divisors));
    
    return divisorsMap[num];
}

void Euler23()
{
    // largest number that can not be written as a sum of two abundant numbers
    int max = 28123;
    
    // create the divisors map
    unordered_map<int, unordered_set<int>> mapOfDivisors;
    for(int i=2; i<max; i++)
    {
        getDivisors(i, mapOfDivisors);
    }
    
    // create array of sum of divisors
    vector<int> sumOfDivisors;
    sumOfDivisors.push_back(0);
    sumOfDivisors.push_back(1);
    
    // array of bools which indicate if this is an abundant number
    vector<bool> isAbundantNumber;
    isAbundantNumber.push_back(false);
    isAbundantNumber.push_back(false);
    
    // array of abundant numbers below max
    vector<int> abundantNumbers;
    
    for(int i=2; i<max; i++)
    {
        unordered_set<int> &divisors = mapOfDivisors[i];
        
        int sum = 0;
        unordered_set<int>::const_iterator divItr = divisors.begin();
        for(; divItr != divisors.end(); ++divItr)
        {
            sum += *divItr;
        }
        
        sumOfDivisors.push_back(sum);
        
        if(sum > i)
        {
            isAbundantNumber.push_back(true);
            abundantNumbers.push_back(i);
        }
        else
            isAbundantNumber.push_back(false);
    }
    
    // using the knowledge that 12 is the smallest abundant number
    int smallestAbundant = 12;
    
    // using the knowledge that 24 is smallest number which can be expressed as a sum of two abundant numbers
    int sumOfnumbers = 23*24/2;     // (adding 1 to 23)
    
    for(int i=25; i<max; i++)
    {
        vector<int>::const_iterator itr = abundantNumbers.begin();
        for(; itr != abundantNumbers.end(); ++itr)
        {
            int secondNumber = i - (*itr);
            
            if(secondNumber < smallestAbundant)
            {
                // this number can not be expressed
                sumOfnumbers += i;
                break;
            }
            
            // if second number is also an abundant number then break
            if(isAbundantNumber[secondNumber])
                break;
        }
    }
    
    cout<<sumOfnumbers<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/

void Euler24()
{
    int maxNum = 10;
    
    // create the array of factorials
    vector<long> factorials;
    factorials.push_back(1);
    for(int i=2; i<maxNum+1; i++)
    {
        factorials.push_back( i*factorials[i-2] );
    }
    
    // this is the lexographic permutation we need to find
    long permutationToFind = 999999;
    
    // these are the digits we are trying to arrange
    vector<int> digits = {0,1,2,3,4,5,6,7,8,9};
    
    // this is the final arrangement of digits
    vector<int> result;
    
    // iterate over the factorials array and find it's location
    for(int i=maxNum-1; i>0; i--)
    {
        long remain = permutationToFind % factorials[i-1];
        long idx= permutationToFind / factorials[i-1];
        
        result.push_back(digits[idx]);
        digits.erase(digits.begin() + idx);
        
        if(remain == 0)
            break;
        
        permutationToFind = remain;
    }
    
    // remove all remaining digits and push them in result
    result.insert(result.end(), digits.begin(), digits.end());
    
    for(int i=0; i<result.size(); i++)
        cout<<result[i];
    
    cout<<endl;
}

/*************************************************************************************************
 *************************************************************************************************/


/*************************************************************************************************
 *************************************************************************************************/


/*************************************************************************************************
 *************************************************************************************************/

void Euler67()
{
    FILE* fp = NULL;
    fp = fopen("triangle67.txt", "r");
    
    int nLines = 100;
    int **nTriangles = new int*[nLines];
    for(int i=0; i<nLines; i++)
    {
        nTriangles[i] = new int[i+1];
        for(int j=0; j<i+1; j++)
            fscanf(fp, "%d", &nTriangles[i][j]);
    }
    
    // this will be a recursive solution with DP
    // allocate memory for holding results
    int **nTriangleSum = new int*[nLines];
    for(int i=0; i<nLines; i++)
    {
        nTriangleSum[i] = new int[i+1];
        for(int j=0; j<i+1; j++)
            nTriangleSum[i][j] = 0;
    }
    
    for(int i=0; i<nLines; i++)
        EvaluateTriangleSum(nTriangles, nTriangleSum, nLines-1, i);
    
    cout<<maxArray(nTriangleSum[nLines-1], nLines)<<endl;
    
    for(int i=0; i<nLines; i++)
    {
        delete [] nTriangles[i];
        delete [] nTriangleSum[i];
    }
    
    delete [] nTriangles;
    delete [] nTriangleSum;
}

/*************************************************************************************************
 *************************************************************************************************/

int main(int argc, const char * argv[]) {
    
    Euler24();
    
    return 0;
}
