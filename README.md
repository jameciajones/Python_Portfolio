# Python_Portfolio
This is the portfolio for python code.

## Analyzing Patient Data
In this analysis we looked at inflammation data for multiple patients.

```python
import numpy
```


```python
numpy.loadtxt(fname = 'inflammation-01.csv', delimiter = ",")
```




    array([[0., 0., 1., ..., 3., 0., 0.],
           [0., 1., 2., ..., 1., 0., 1.],
           [0., 1., 1., ..., 2., 1., 1.],
           ...,
           [0., 1., 1., ..., 1., 1., 1.],
           [0., 0., 0., ..., 0., 2., 0.],
           [0., 0., 1., ..., 1., 1., 0.]])




```python
data = numpy.loadtxt(fname = 'inflammation-01.csv', delimiter = ",")
```


```python
print(data)
```

    [[0. 0. 1. ... 3. 0. 0.]
     [0. 1. 2. ... 1. 0. 1.]
     [0. 1. 1. ... 2. 1. 1.]
     ...
     [0. 1. 1. ... 1. 1. 1.]
     [0. 0. 0. ... 0. 2. 0.]
     [0. 0. 1. ... 1. 1. 0.]]



```python
print(data.shape)
```

    (60, 40)



```python
print('first value in data:' , data [0,0])
```

    first value in data: 0.0



```python
print('middle value in data:' , data[29,19])
```

    middle value in data: 16.0



```python
print(data[0:4, 0:10])
```

    [[0. 0. 1. 3. 1. 2. 4. 7. 8. 3.]
     [0. 1. 2. 1. 2. 1. 3. 2. 2. 6.]
     [0. 1. 1. 3. 3. 2. 6. 2. 5. 9.]
     [0. 0. 2. 0. 4. 2. 2. 1. 6. 7.]]



```python
print(data[0:5, 0:10])
```

    [[0. 0. 1. 3. 1. 2. 4. 7. 8. 3.]
     [0. 1. 2. 1. 2. 1. 3. 2. 2. 6.]
     [0. 1. 1. 3. 3. 2. 6. 2. 5. 9.]
     [0. 0. 2. 0. 4. 2. 2. 1. 6. 7.]
     [0. 1. 1. 3. 3. 1. 3. 5. 2. 4.]]



```python
small = data[:3, 36:]
```


```python
print('small is:')
```

    small is:



```python
print(small)
```

    [[2. 3. 0. 0.]
     [1. 1. 0. 1.]
     [2. 2. 1. 1.]]



```python
# Let us a numpy function
print(numpy.mean(data))
```

    6.14875



```python
maxval, minval, stdvat = numpy.amax(data), numpy.amin(data), numpy.std(data)


```


```python
print(maxval)
print(minval)
print(stdval)
```

    20.0
    0.0



    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-30-566e6d04f418> in <module>
          1 print(maxval)
          2 print(minval)
    ----> 3 print(stdval)
    

    NameError: name 'stdval' is not defined



```python
maxval = numpy.amax(data)
minval = numpy.amin(data)
stdval = numpy.std(data)
```


```python
print(maxval)
print(minval)
print(stdval)
```

    20.0
    0.0
    4.613833197118566



```python
print('maximum inflammation:' , maxval)
print('minimum inflammation:' , minval)
print('standard deviation:' , stdval)
```


```python
# Sometimes we want to look at variation in statistical values, such as maximum inflammation per patient, 
# or average from day one

patient_0 = data[0, :] # 0 on the first axis (rows), everything on the second (columns)

print('maximum inflammation for patient 0:', numpy.amax(patient_0))
```


```python
print('maximum inflammation for patient:2' , numpy.amax(data[2, :]))
```


```python
print(numpy.mean(data, axis =0))
```


```python
print(numpy.mean(data, axis = 0).shape)
```


```python
print(numpy.mean(data, axis = 1))
```


## Storing Values in Lists


```python
odds = [1, 3, 5, 7]
print('odds are:', odds)
```

    odds are: [1, 3, 5, 7]



```python
print('first element:', odds [0])
print('last element:' , odds [3])
print('"-1" element:', odds [-1])
```

    first element: 1
    last element: 7
    "-1" element: 7



```python
names = ['Curie' ,'Darwing', 'Turing'] # Typo in Darwin's name

print('names is originally:' , names)

name[1] = 'Darwin' # Correct the name

print('final value of names:' , names)
```

    names is originally: ['Curie', 'Darwing', 'Turing']
    final value of names: ['Curie', 'Darwing', 'Turing']



```python
#name = 'Darwin'
#name[0] = 'd'
```


```python
odds.append(11)
print('odds after adding a value:' , odds)
```

    odds after adding a value: [1, 3, 5, 7, 11]



```python
removed_element = odds.pop(0)
print('odds after removing the first element:' , odds)
print('removed_element:' , removed_element)
```

    odds after removing the first element: [3, 5, 7, 11]
    removed_element: 1



```python
odds.reverse()
print('odds after reversing:' , odds)
```

    odds after reversing: [11, 7, 5, 3]



```python
odds = [3,5,7]
primes = odds
primes.append(2)
print('primes;' , primes)
print('odds:' , odds)
```

    primes; [3, 5, 7, 2]
    odds: [3, 5, 7, 2]



```python
odds = [3,5,7]
primes = list(odds)
primes.append(2)
print('primes;' , primes)
print('odds:' , odds)
```

    primes; [3, 5, 7, 2]
    odds: [3, 5, 7]



```python
binomial_name = "Drosophila melangogaster"
group = binomial_name[0:10]
print('group:' , group)

species = binomial_name[11:23]
print('species:' , species)

chromosomes = ['X' , 'Y' , '2' , '3', '4']
autosomes = chromosomes[2:5]
print('autosomes:' , autosomes)

last = chromosomes [-1]
print('last:', last)
```

    group: Drosophila
    species: melangogaste
    autosomes: ['2', '3', '4']
    last: 4



```python
date = 'Monday 4 January 2023'
day = date[0:6]
print('Using 0 to begin range:', day)
day = date[:6]
print('Omitting beginning index:', day)
```

    Using 0 to begin range: Monday
    Omitting beginning index: Monday



```python
months = ['jan', 'feb' , 'mar' , 'apr' ,'may' , 'jun', 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec']
sond = months

print('With known last postion:' , sond)

sond = months[8:len(months)]
print('Using len() to get last entry:' , sond)

sond = months[8:]
print('Omitting ending index:' , sond)


```

    With known last postion: ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
    Using len() to get last entry: ['sep', 'oct', 'nov', 'dec']
    Omitting ending index: ['sep', 'oct', 'nov', 'dec']

## Using Loops

```python
odds = [1,3,5,7]
```


```python
print(odds[0])
print(odds[1])
print(odds[2])
print(odds[3])
```

    1
    3
    5
    7



```python
odds = [1,3,5]
print(odds[0])
print(odds[1])
print(odds[2])
print(odds[3])
```

    1
    3
    5



    ---------------------------------------------------------------------------

    IndexError                                Traceback (most recent call last)

    <ipython-input-5-01ba67d8a9e5> in <module>
          3 print(odds[1])
          4 print(odds[2])
    ----> 5 print(odds[3])
    

    IndexError: list index out of range



```python
odds = [1,3,5,7,9,11.13,15,17,19]

for num in odds:
    print(num)
```

    1
    3
    5
    7
    9
    11.13
    15
    17
    19



```python
lenght = 0
names = ['Curie,' 'Darwin', 'Turing']
for value in names:
    lenght = lenght + 1
print('There are', lenght, 'names in the list.')

```

    There are 2 names in the list.



```python
name = "Rosalind"
for name in ['Curie', 'Darwin' , 'Turing']:
    print(name)
print('after the loop, name is' , name)
```

    Curie
    Darwin
    Turing
    after the loop, name is Turing



```python
print(len([0,1,2,3]))
```

    4



```python
name = ['Curie,' 'Darwin', 'Turing']

print(len(name))
```

    2

## Python Fundamentals

```python
# Any python interpreter can be used as a calculator:
3 + 5 + 4

```




    12




```python
# Lets save a value to a variable
weight_kg = 60 
```


```python
print(weight_kg)
```

    60



```python
# Weight0 = valud
# 0weight = invalid
# weight and Weight are different 
```


```python
# Types of data
# There are three common types of data
# Integer numbers
# floating point numbers
# Strings

```


```python
# Floating point number
weight_kg = 60.3

```


```python
# String comprised of letters
patient_name = "Jon Smith"
```


```python
# String comprised of number
patient_id = '001'
```


```python
# Use variables in python 

weight_lb = 2.2 * weight_kg

print(weight_lb)
```

    132.0



```python
# Lets add a prefix to our patient id 

patient_id = 'inflam-' + patient_id

print(patient_id)
```

    inflam-001



```python
# Lets combine print statements 

print(patient_id, 'weight in kilogram:', weight_kg)
```

    inflam-001 weight in kilogram: 60



```python
# we can call a function inside another function

print(type(60.3))

print(type(patient_id))
```

    <class 'float'>
    <class 'str'>



```python
# We can also do calculation inside the print function

print('weight in lbs:', 2.2 * weight_kg)
```

    weight in lbs: 132.0



```python
print(weight_kg)
```

    60



```python
weight_kg = 65.0
print('weight in kilograms is now:', weight_kg)
```

    weight in kilograms is now: 65.0

## Using Multiple Files
```python
import glob

```


```python
print(glob.glob('inflammation*.csv'))
```

    ['inflammation-05.csv', 'inflammation-12.csv', 'inflammation-04.csv', 'inflammation-08.csv', 'inflammation-10.csv', 'inflammation-06.csv', 'inflammation-09.csv', 'inflammation-01.csv', 'inflammation-07.csv', 'inflammation-11.csv', 'inflammation-03.csv', 'inflammation-02.csv']



```python
import glob
import numpy
import matplotlib.pyplot

filenames = sorted(glob.glob('inflammation*.scv'))
filesnames = filenames[0:3]

for filename in filenames:
    print(filename)
    
    data = numpy.loadtxt(fname=filename, delimiter = ',')
    
    fig = matplotlib.pyplot.figure(figsize = (10.0 , 3.0))
    
    axes1 = fig.add_subplot(1,3,1)
    axes2 = fig.add_subplot(1,3,2)
    axes3 = fig.add_subplot(1,3,3)
    
    axes1.set_ylabel('average')
    axes1.plot(numpy.mean(data,axis = 0))
  
    axes2.set_ylabel('max')
    axes2.plot(numpy.amax(data, axis = 0))
    
    axes3.set_ylabel('min')
    axes3.plot(numpy.amin(data,axis = 0))
    
    fig.tight_layout()
    matplotlib.pyplot.show()
```

## Making Choices

```python
num = 37
if num > 100:
    print('greater')
else:
    print('not greater')
print('done')
    
```

    not greater
    done



```python
num = 53
print('before conditional...')
if num > 100:
   print(num, 'is greater than 100')
print('...after conditional')
```

    before conditional...
    ...after conditional


```python
nuum = -3

if num > 0:
    print(num,'is postive')
elif num == 0:
    print(num, 'is zero')
else: 
   
## Functions
