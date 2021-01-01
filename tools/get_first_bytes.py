# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: get_first_bytes.py
# time: 2020/10/24
from argparse import ArgumentParser
from ftplib import FTP
import os


class LimitedFile:
    def __init__(self, filename, size=512):
        self.filename = filename
        self.size = size
        self.count = 0
        if os.path.isfile(filename):
            os.remove(filename)

    def write(self, data):
        data = data[:self.size]
        with open(self.filename, 'ab') as f:
            c = f.write(data)
        self.count += c
        if self.count >= self.size:
            raise ValueError("Reached limitation")


def main():
    parser = ArgumentParser(
        description='get the file\'s first bytes from a ftp server'
    )
    parser.add_argument('<url>', help='ftp link')
    parser.add_argument('-c', default=512, help='size you want')
    args = vars(parser.parse_args())
    host, *directory, filename = args['<url>'].split('/')
    ftp = FTP(host)
    ftp.login()
    f = LimitedFile(filename, args['c'])
    try:
        ftp.retrbinary(f'RETR /{"/".join(directory)}/{filename}', f.write, 512)
    except ValueError:
        ftp.abort()
    ftp.close()
    print(f'Get {f.count} bytes from {args["<url>"]}')


if __name__ == '__main__':
    main()
