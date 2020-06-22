class Gene:

    def __init__(self, x_cor, y_cor):
        self.x_cor = x_cor
        self.y_cor = y_cor

    def equals(self, point):
        is_same = False
        if self.x_cor == point.x_cor and self.y_cor == point.y_cor:
            is_same = True
        return is_same

    def print(self, message):
        print(message, " (", self.x_cor, ", ", self.y_cor, ")")



