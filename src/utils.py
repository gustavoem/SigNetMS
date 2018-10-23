def safe_power (a, b):
    try:
        x = a ** b
    except OverflowError:
        x = float ("inf")
    return x
