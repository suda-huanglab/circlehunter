# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: test.py
# time: 2020/11/19
from multiprocessing import Process
from queue import Queue, Empty
from functools import partial
import time


def worker(queue):
    print('worker init')
    while True:
        try:
            queue.get_nowait()
        except Empty:
            print('queue is empty')
            continue
        print('got item')
        time.sleep(1)


def run():
    queue = Queue(maxsize=10)
    target = partial(time.sleep, 1)
    workers = [Process(target=target) for _ in range(2)]
    for process in workers:
        process.start()
    while True:
        if queue.full():
            # print('queue is full')
            continue
        queue.put(1)


if __name__ == '__main__':
    run()
