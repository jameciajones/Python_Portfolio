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
## Functions
```python
fahrenheit_val = 99
celsius_val = ((fahrenheit_val - 32) * (5/9))
             
print(celsius_val)
```

    37.22222222222222



```python
fahrenheit_val = 43
celsius_val2 = ((fahrenheit_val - 32) * (5/9))
             
print(celsius_val2)
```

    6.111111111111112



```python
def explicit_fahr_celcius(temp):
   # Assign the converted value to a variable
   converted = ((temp - 32) * (5/9))
   # Return the values of the new variable
   return converted


```


```python
def fahr_to_celsius(temp):
    # Return converted value more ffeciently using the return function without creating
    # a new variable. This code does the same thing as the previous function but it is more
    # explicit in explaining how the return command works.
    return ((temp - 32) * (5/9))
```


```python
fahr_to_celsius(32)
```




    0.0




```python
print('Freezing point of water:', fahr_to_celsius(32), 'C')
print('Boiling point of water:', fahr_to_celsius(212))
```

    Freezing point of water: 0.0 C
    Boiling point of water: 100.0



```python
def celsius_to_kelvin(temp_c):
    return temp_c + 273.15
    
print('freezing point of water in Kelvin:', celsius_to_kelvin(0.))
```

    freezing point of water in Kelvin: 273.15



```python
def fahr_to_kelvin(temp_f):
    temp_c = fahr_to_celsius(temp_f)
    temp_k = celsius_to_kelvin(temp_c)
    return temp_k

print('boiling point of water in Kelvin:', fahr_to_kelvin(212.0))
```

    boiling point of water in Kelvin: 373.15



```python
print('Again, temperature in Kelving was:', temp_k)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-37-be27d2dc3254> in <module>
    ----> 1 print('Again, temperature in Kelving was:', temp_k)
    

    NameError: name 'temp_k' is not defined



```python
temp_kelving = fahr_to_kelvin(212.0)
print('Temperature in Kelvin was:', temp_kelving)
```

    Temperature in Kelvin was: 373.15



```python
temp_kelvin
```




    373.15




```python
def print_temperature():
    print('Temperature in Fahrenheit was:', temp_fahr)
    print('Temperature in Kelvin was:', temp_kelvin)

temp_fahr = 212.0
temp_kelvin = fahr_to_kelvin(temp_fahr)

print_temperature()
```

    Temperature in Fahrenheit was: 212.0
    Temperature in Kelvin was: 373.15



```python

```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/87771406-c8d6-4701-927c-ab07399524fc)

```python
import numpy 
import matplotlib
import matplotlib.pyplot
import glob
```


```python
'freezing point of water in Kelvin:'
def visualize (filename):
    
    data = numpy.loadtxt(fname = filename, delimiter = ",")
    
    fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))
    
    axes1 = fig.add_subplot(1, 3, 1)
    axes2 = fig.add_subplot(1, 3, 2)
    axes3 = fig.add_subplot(1, 3, 3)
    
    axes1.set_ylabel('average')
    axes1.plot(numpy.mean(data,axis=0))
    
    
    axes2.set_ylabel('max')
    axes2.plot(numpy.amax(data, axis = 0))
    
    axes3.set_ylabel('min')
    axes3.plot(numpy.amin(data, axis = 0))
    
    fig.tight_layout()
    matplotlib.pyplot.show()
```


```python
def detect_problems(filename):
    
    data = numpy.loadtxt(fname = filename, delimiter = ",")
    
    if numpy.amax(data, axis = 0)[0] == 0 and numpy.amax(data, axis =0)[20] == 20:
        print("Suspicious looking maxima!")
    elif numpy.sum(numpy.amin(data, axis=0)) == 0:
        print('Minima add up to zero')
    else:
        print('Seems ok!')
```


```python
filenames = sorted(glob.glob('inflammation*.csv'))

for filename in filenames[:3]:
    print(filename)
    visualize(filename)
    detect_problems(filename)
```

    inflammation-01.csv



![png](output_3_1.png)


    Suspicious looking maxima!
    inflammation-02.csv



![png](output_3_3.png)


    Suspicious looking maxima!
    inflammation-03.csv



![png](output_3_5.png)


    Minima add up to zero



```python
def offset_mean(data, target_mean_value):
    return(data - numpy.mean(data)+ target_mean_value)
```


```python
z = numpy.zeros((2,2))
print(offset_mean(z,3))
```

    [[3. 3.]
     [3. 3.]]



```python
data = numpy.loadtxt(fname = 'inflammation-01.csv', delimiter = ',')

print(offset_mean(data,0))
```

    [[-6.14875 -6.14875 -5.14875 ... -3.14875 -6.14875 -6.14875]
     [-6.14875 -5.14875 -4.14875 ... -5.14875 -6.14875 -5.14875]
     [-6.14875 -5.14875 -5.14875 ... -4.14875 -5.14875 -5.14875]
     ...
     [-6.14875 -5.14875 -5.14875 ... -5.14875 -5.14875 -5.14875]
     [-6.14875 -6.14875 -6.14875 ... -6.14875 -4.14875 -6.14875]
     [-6.14875 -6.14875 -5.14875 ... -5.14875 -5.14875 -6.14875]]



```python
print('original min, mean and max are:', numpy.amin(data), numpy.mean(data), numpy.amax(data))
offset_data = offset_mean(data, 0)
print('min, mean, and max of offset data are:' ,
      numpy.amin(offset_data),
      numpy.mean(offset_data),
      numpy.amax(offset_data))
```

    original min, mean and max are: 0.0 6.14875 20.0
    min, mean, and max of offset data are: -6.14875 2.842170943040401e-16 13.85125



```python
print('std dev before and after:', numpy.std(data), numpy.std(offset_data))
```

    std dev before and after: 4.613833197118566 4.613833197118566



```python
print('difference in standard deviation before and after:',
     numpy.std(data) - numpy.std(offset_data))
```

    difference in standard deviation before and after: 0.0



```python
# offset_mean(data, target_mean_value):
# return a new array containing the orginal data with its mean offset to match the desired value.

def offset_mean(data, target_mean_value):
    return (data - numpy.mean(data)) + target_mean_value
```


```python
def offset_mean(data, target_mean_value):
    """"Return a new arrary containg the original datra with its mean offset to match th desired value"""
    return(data - numpy.mean(data)) + target_mean_value
    
```


```python
help(offset_mean)
```

    Help on function offset_mean in module __main__:
    
    offset_mean(data, target_mean_value)
        "Return a new arrary containg the original datra with its mean offset to match th desired value
    



```python
def offset_mean(data, target_mean_value):
    """"Return a new arrary containing the orginal data with its mean offset to match the desired value.
    
    Examples
    --------
    
    >>> Offset_mean([1,2,3],0)
    array([1., 0., 1])
    
    """
    
    return (data - numpy.mean(data)) + target_mean_value
```


```python
help(offset_mean)
```

    Help on function offset_mean in module __main__:
    
    offset_mean(data, target_mean_value)
        "Return a new arrary containing the orginal data with its mean offset to match the desired value.
        
        Examples
        --------
        
        >>> Offset_mean([1,2,3],0)
        array([1., 0., 1])
    



```python
numpy.loadtxt('inflammation-01.csv', delimiter = ',')

```




    array([[0., 0., 1., ..., 3., 0., 0.],
           [0., 1., 2., ..., 1., 0., 1.],
           [0., 1., 1., ..., 2., 1., 1.],
           ...,
           [0., 1., 1., ..., 1., 1., 1.],
           [0., 0., 0., ..., 0., 2., 0.],
           [0., 0., 1., ..., 1., 1., 0.]])




```python
def offset_mean(data,target_mean_value = 0.0):
    """"Return a new arrary containing the orginal data with its mean offset to match the desired value,
    (0 by default).
    
    Examples
    ---------
    
    >>> offset_mean([1,2,3], 0)
    array ([-1.,0.,1.])
    
    """
    
    return(data - numpy.mean(data) + target_mean_value)
```


```python
test_data = numpy.zeros((2,2))
print(offset_mean(test_data, 3))
```

    [[3. 3.]
     [3. 3.]]



```python
print(offset_mean(test_data))
```

    [[0. 0.]
     [0. 0.]]



```python
def display (a=1, b=2, c=3):
    print('a:', a, 'b:', b, 'c:', c)
    
print('no parameters:')
display()
print('one parameter:')
display(55)
print('two parameters:')
display(55, 66)
```

    no parameters:
    a: 1 b: 2 c: 3
    one parameter:
    a: 55 b: 2 c: 3
    two parameters:
    a: 55 b: 66 c: 3



```python
print('only setting the value of c')
display(c = 77)
```

    only setting the value of c
    a: 1 b: 2 c: 77



```python
help(numpy.loadtxt)
```

    Help on function loadtxt in module numpy:
    
    loadtxt(fname, dtype=<class 'float'>, comments='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0, encoding='bytes', max_rows=None)
        Load data from a text file.
        
        Each row in the text file must have the same number of values.
        
        Parameters
        ----------
        fname : file, str, or pathlib.Path
            File, filename, or generator to read.  If the filename extension is
            ``.gz`` or ``.bz2``, the file is first decompressed. Note that
            generators should return byte strings for Python 3k.
        dtype : data-type, optional
            Data-type of the resulting array; default: float.  If this is a
            structured data-type, the resulting array will be 1-dimensional, and
            each row will be interpreted as an element of the array.  In this
            case, the number of columns used must match the number of fields in
            the data-type.
        comments : str or sequence of str, optional
            The characters or list of characters used to indicate the start of a
            comment. None implies no comments. For backwards compatibility, byte
            strings will be decoded as 'latin1'. The default is '#'.
        delimiter : str, optional
            The string used to separate values. For backwards compatibility, byte
            strings will be decoded as 'latin1'. The default is whitespace.
        converters : dict, optional
            A dictionary mapping column number to a function that will parse the
            column string into the desired value.  E.g., if column 0 is a date
            string: ``converters = {0: datestr2num}``.  Converters can also be
            used to provide a default value for missing data (but see also
            `genfromtxt`): ``converters = {3: lambda s: float(s.strip() or 0)}``.
            Default: None.
        skiprows : int, optional
            Skip the first `skiprows` lines, including comments; default: 0.
        usecols : int or sequence, optional
            Which columns to read, with 0 being the first. For example,
            ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns.
            The default, None, results in all columns being read.
        
            .. versionchanged:: 1.11.0
                When a single column has to be read it is possible to use
                an integer instead of a tuple. E.g ``usecols = 3`` reads the
                fourth column the same way as ``usecols = (3,)`` would.
        unpack : bool, optional
            If True, the returned array is transposed, so that arguments may be
            unpacked using ``x, y, z = loadtxt(...)``.  When used with a structured
            data-type, arrays are returned for each field.  Default is False.
        ndmin : int, optional
            The returned array will have at least `ndmin` dimensions.
            Otherwise mono-dimensional axes will be squeezed.
            Legal values: 0 (default), 1 or 2.
        
            .. versionadded:: 1.6.0
        encoding : str, optional
            Encoding used to decode the inputfile. Does not apply to input streams.
            The special value 'bytes' enables backward compatibility workarounds
            that ensures you receive byte arrays as results if possible and passes
            'latin1' encoded strings to converters. Override this value to receive
            unicode arrays and pass strings as input to converters.  If set to None
            the system default is used. The default value is 'bytes'.
        
            .. versionadded:: 1.14.0
        max_rows : int, optional
            Read `max_rows` lines of content after `skiprows` lines. The default
            is to read all the lines.
        
            .. versionadded:: 1.16.0
        
        Returns
        -------
        out : ndarray
            Data read from the text file.
        
        See Also
        --------
        load, fromstring, fromregex
        genfromtxt : Load data with missing values handled as specified.
        scipy.io.loadmat : reads MATLAB data files
        
        Notes
        -----
        This function aims to be a fast reader for simply formatted files.  The
        `genfromtxt` function provides more sophisticated handling of, e.g.,
        lines with missing values.
        
        .. versionadded:: 1.10.0
        
        The strings produced by the Python float.hex method can be used as
        input for floats.
        
        Examples
        --------
        >>> from io import StringIO   # StringIO behaves like a file object
        >>> c = StringIO(u"0 1\n2 3")
        >>> np.loadtxt(c)
        array([[0., 1.],
               [2., 3.]])
        
        >>> d = StringIO(u"M 21 72\nF 35 58")
        >>> np.loadtxt(d, dtype={'names': ('gender', 'age', 'weight'),
        ...                      'formats': ('S1', 'i4', 'f4')})
        array([(b'M', 21, 72.), (b'F', 35, 58.)],
              dtype=[('gender', 'S1'), ('age', '<i4'), ('weight', '<f4')])
        
        >>> c = StringIO(u"1,0,2\n3,0,4")
        >>> x, y = np.loadtxt(c, delimiter=',', usecols=(0, 2), unpack=True)
        >>> x
        array([1., 3.])
        >>> y
        array([2., 4.])
    



```python
numpy.loadtxt('inflammation-01.csv', delimiter = ',')
```




    array([[0., 0., 1., ..., 3., 0., 0.],
           [0., 1., 2., ..., 1., 0., 1.],
           [0., 1., 1., ..., 2., 1., 1.],
           ...,
           [0., 1., 1., ..., 1., 1., 1.],
           [0., 0., 0., ..., 0., 2., 0.],
           [0., 0., 1., ..., 1., 1., 0.]])




```python

```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/c4830db8-34de-49e1-95be-042e9b3edeb0)
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
    
```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/24f5d724-dd40-40f5-9007-0836b05b916d)
## Errors
```python
# This code has an intentional error. You can type it directly or use it for reference 
# to understand the error message below.

def favorite_ice_cream():
    ice_creams = [
        'chacolate'
        'vanilla'
        'strawberry'
    ]
    
    print(ice_creams[3])
    
    favorite_ice_cream()
```


```python
def some_function():
    msg = 'hello,world!'
    print(msg)
    return msg
```


```python
print(a)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-5-bca0e2660b9f> in <module>
    ----> 1 print(a)
    

    NameError: name 'a' is not defined



```python
print('hello')
```

    hello



```python
count = 0

for number in range(10):
    count = count + number
print('The count is:', count)
```

    The count is: 45



```python
letters = ['a', 'b', 'c']

print('Letter #1 is', letters[0])
print('Letter #2 is', letters [1])
print('Letter #3 is', letters[2])
print('Letter #4 is', letters[2])
```

    Letter #1 is a
    Letter #2 is b
    Letter #3 is c
    Letter #4 is c



```python
file_handle = open('myfile.txt', 'w')
```


```python

```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/403c71e3-fa00-49da-886d-347b22371c85)
## Defensive Programming
```python
numbers = [1.5, 2.3, 0.7, 0.001, 4.4]
total = 0.0
for num in numbers:
    assert num > 0.0, 'Data should only contain positive valuse'
    total += num
print('total is:', total)
```

    total is: 8.901



```python
def normalize_rectangle(rect):
    """Normalizes a rectangle so that it is at the orgin and 1.0 units long on its
    logest axis. input should be of the format of (x0,y0, x1,y1).
    (x0,y0) and (x1,y1) define the lower left and upper right corners of the 
    rectangle respectively."""
    assert len(rect) == 4, 'Rectangles must contain 4 coordniates'
    x0,y0,x1,y1 = rect
    assert x0 < x1, 'Invalid X coordinates'
    assert y0 < y1, 'Invalid Y coordinates'
    
    dx = x1 - x0
    dy = y1 - y0
    if dx > dy:
        scaled = dx / dy
        upper_x, upper_y = 1.0, scaled
    else:
        scaled = dx / dy
        upper_x, upper_y = scaled, 1.0
        
    assert 0 < upper_x <= 1.0, 'Calculated upper x coordinate invalid'
    assert 0 < upper_y <= 1.0, 'Calculated upper y coordinate invalid'
    
    return (0, 0, upper_x, upper_y)
```


```python
print(normalize_rectangle((0.0, 0.0, 1.0, 5.0)) )
```

    (0, 0, 0.2, 1.0)



```python

```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/9de80f55-394d-4af9-8339-6cc47f50b261)
## Command Line Programs
```python
import numpy
data = numpy.loadtxt(fname= 'inflammation-01.csv' , delimiter = ',')
```


```python
import matplotlib.pyplot
image = matplotlib.pyplot.imshow(data)
matplotlib.pyplot.show()
```


![png](output_1_0.png)



```python
# Average inflammation over time

ave_inflammation = numpy.mean(data, axis = 0)
ave_plot = matplotlib.pyplot.plot(ave_inflammation)
matplotlib.pyplot.show()
```


![png](output_2_0.png)



```python
max_plot = matplotlib.pyplot.plot(numpy.amax(data, axis = 0))
matplotlib.pyplot.show()
```


![png](output_3_0.png)



```python
min_plot = matplotlib.pyplot.plot(numpy.amin(data, axis = 0))
matplotlib.pyplot.show()
```


![png](output_4_0.png)



```python
fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))

axes1 = fig.add_subplot(1, 3, 1)
axes2 = fig.add_subplot(1, 3, 2)
axes3 = fig.add_subplot(1, 3, 3)

axes1.set_ylabel('average')
axes1.plot(numpy.mean(data, axis =0))

axes2.set_ylabel('max')
axes2.plot(numpy.amax(data, axis = 0))

axes3.set_ylabel('min')
axes3.plot(numpy.amin(data, axis =0))

fig.tight_layout()

matplotlib.pyplot.savefig('inflammation.png')
matplotlib.pyplot.show()
```


![png](output_5_0.png)



```python

```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/08e1600b-2da0-41f9-9873-90fb134738fb)
## Transcribing DNA to RNA
```python
# Prompt the user to enter the input fasta file name 

input_file_name = input("Enter the name of the input fasta file: ")
```

    Enter the name of the input fasta file:  SUMO.txt



```python
# Open the input fasta file and read the DNA sequence 

with open(input_file_name, "r") as input_file:
    dna_sequence = ""
    for line in input_file:
        if line.startswith(">"):
            continue
        dna_sequence += line.strip()
```


```python
# Transcribe the DNA to RNA
rna_sequence = ""
for nucleotide in dna_sequence:
    if nucleotide == "T":
        rna_sequence += "U"
    else:
        rna_sequence += nucleotide
```


```python
# Prompt the yser to enter the output file name

output_file_name = input("Enter the name of the output file: ")
```

    Enter the name of the output file:  Ubiquitin_RNA.txt



```python
# Save the RNA sequence to a text file

with open(output_file_name, "w") as output_file:
    output_file.write(rna_sequence)
    print("The RNA sequence has been saved to {output_file_name}")
```

    The RNA sequence has been saved to {output_file_name}



```python
print(rna_sequence)
```

    AUGUCUGACGAAAAGAAGGGAGGUGAGACCGAGCACAUCAACCUGAAGGUCCUCGGCCAGGACAACGCCGUCGUCCAGUUCAAGAUCAAGAAGCACACACCCUUGAGGAAGCUGAUGAACGCCUACUGCGACCGUGCCGGACUCUCCAUGCAGGUGGUGCGCUUCCGUUUCGACGGACAGCCCAUCAACGAGAACGACACUCCGACCUCGCUGGAGAUGGAGGAGGGCGACACCAUCGAGGUUUACCAGCAGCAGACUGGUGGCGCUCCAUAAAUGUCUGACGAAAAGAAGGGAGGUGAGACCGAGCACAUCAACCUGAAGGUCCUCGGCCAGGACAACGCCGUCGUCCAGUUCAAGAUCAAGAAGCACACACCCUUGAGGAAGCUGAUGAACGCCUACUGCGACCGUGCCGGACUCUCCAUGCAGGUGGUGCGCUUCCGUUUCGACGGACAGCCCAUCAACGAGAACGACACUCCGACCUCGCUGGAGAUGGAGGAGGGCGACACCAUCGAGGUUUACCAGCAGCAGACUGGUGGCGCUCCAUAA



```python


```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/4b7e733b-f6db-4c45-bd66-9d3191937205)
## Translating RNA into Protein
```python
# Prompt the user to enter the input RNA file name

input_file_name = input("Enter the name of the input RNA file:")
```

    Enter the name of the input RNA file: Ubiquitin_RNA.txt



```python
# Open the input RNA file and read the RNA sequence

with open(input_file_name, "r") as input_file:
    rna_sequence = input_file.read().strip()
```


```python
# Define the codon table 

codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UCA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    
}
```


```python
# Translate RNA to protein

protein_sequence = " "
for i in range(0,len(rna_sequence), 3):
    if len(codon) == 3:
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        protein_sequence += amino_acid
   
```


```python
# Prompt the user to enter the output file name

output_file_name = input("Enter the name of the output file: ")
```

    Enter the name of the output file:  Ubiquitin_Protein.txt



```python
# Save the protein sequence to a text file

with open(output_file_name, "w") as output_file:
    output_file.write(protein_sequence)
    print(f"The protein sequence has been saved to {output_file_name}")
    
```

    The protein sequence has been saved to Ubiquitin_Protein.txt



```python
print(protein_sequence)
```

      MSDEKKGGETEHINLKVLGQDNAVVQFKIKKHTPLRKLMNAYCDRAGLSMQVVRFRFDGQPINENDTPTSLEMEEGDTIEVYQQQTGGAP



```python

```

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/8ba5fb83-86d1-4ab3-9a62-fe8b186b9ee9)
## Using Jupyter Notebooks
```python %matplotlib inline import pandas as pd import matplotlib.pyplot as plt import seaborn as sns sns.set(style="darkgrid") ``` ```python df = pd.read_csv('/home/student/Desktop/classroom/myfiles/notebooks/fortune500.csv') ``` ```python df.head() ``` 
	Year	Rank	Company	Revenue (in millions)	Profit (in millions)
0	1955	1	General Motors	9823.5	806
1	1955	2	Exxon Mobil	5661.4	584.8
2	1955	3	U.S. Steel	3250.4	195.4
3	1955	4	General Electric	2959.1	212.6
4	1955	5	Esmark	2510.8	19.1
```python df.tail() ``` 
	Year	Rank	Company	Revenue (in millions)	Profit (in millions)
25495	2005	496	Wm. Wrigley Jr.	3648.6	493
25496	2005	497	Peabody Energy	3631.6	175.4
25497	2005	498	Wendy's International	3630.4	57.8
25498	2005	499	Kindred Healthcare	3616.6	70.6
25499	2005	500	Cincinnati Financial	3614.0	584
```python df.columns = ['year', 'rank','company', 'revenue','profit'] ``` ```python df.head() ``` 
	year	rank	company	revenue	profit
0	1955	1	General Motors	9823.5	806
1	1955	2	Exxon Mobil	5661.4	584.8
2	1955	3	U.S. Steel	3250.4	195.4
3	1955	4	General Electric	2959.1	212.6
4	1955	5	Esmark	2510.8	19.1
```python len(df) ``` 25500 ```python df.dtypes ``` year int64 rank int64 company object revenue float64 profit object dtype: object ```python non_numeric_profits = df.profit.str.contains('[^0-9.-]') df.loc[non_numeric_profits].head() ``` 
	year	rank	company	revenue	profit
228	1955	229	Norton	135.0	N.A.
290	1955	291	Schlitz Brewing	100.0	N.A.
294	1955	295	Pacific Vegetable Oil	97.9	N.A.
296	1955	297	Liebmann Breweries	96.0	N.A.
352	1955	353	Minneapolis-Moline	77.4	N.A.
```python set(df.profit[non_numeric_profits]) ``` {'N.A.'} ```python len(df.profit[non_numeric_profits]) ``` 369 ```python bin_sizes, _,_ = plt.hist(df.year[non_numeric_profits], bins=range(1955,2006) ) ``` ![png](output_11_0.png) ```python df = df.loc[~non_numeric_profits] df.profit = df.profit.apply(pd.to_numeric) ``` --------------------------------------------------------------------------- NameError Traceback (most recent call last) in ----> 1 df = df.loc[~non_numeric_profits] 2 df.profit = df.profit.apply(pd.to_numeric) NameError: name 'df' is not defined ```python len(df) ``` --------------------------------------------------------------------------- NameError Traceback (most recent call last) in ----> 1 len(df) NameError: name 'df' is not defined ```python df = df.loc[~non_numeric_profits] df.profit = df.profit.apply(pd.to_numeric) ``` --------------------------------------------------------------------------- NameError Traceback (most recent call last) in ----> 1 df = df.loc[~non_numeric_profits] 2 df.profit = df.profit.apply(pd.to_numeric) NameError: name 'df' is not defined ```python len (df) ``` --------------------------------------------------------------------------- NameError Traceback (most recent call last) in ----> 1 len (df) NameError: name 'df' is not defined ```python a ``` 
![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/d2e5d9d7-fd06-452b-b16b-46884de75f15)


    
    
 
   


