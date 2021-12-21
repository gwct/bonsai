import sys, os
import multiprocessing as mp

def square(x):
    created = mp.Process();
    current = mp.current_process();
    print('running:', current.name, current._identity);
    print('created:', created.name, created._identity);
    y = list(range(0, 300000, 2));

    if x in y:
        return "NOPE";

    else:
        return x * x;


nums = list(range(0, 300000));

# for num in nums:
#     print(square(num));

pool = mp.Pool(processes=4);

with open("test.txt", "a") as outfile:
    #result = pool.imap(square, nums);
    for result in pool.imap(square, nums):
        print(result);
        #if result is not None:
        outfile.write(str(result) + "\n");
        outfile.flush();
        #print(result);
    #print(len(result));


pool.close();