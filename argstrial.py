import argparse


def check_time(value):
    ivalue = int(value)
    if ivalue > 90 or ivalue < 45:
        raise argparse.ArgumentTypeError("%s is an invalid time value" % value)
    return ivalue


if __name__  == '__main__' :

    #parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(...)
    parser.add_argument('foo', type=check_time)
    args = parser.parse_args()


