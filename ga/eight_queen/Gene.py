class Gene:

    def __init__(self, x_cor, y_cor):
        if not isinstance(x_cor, int) or not isinstance(y_cor, int):
            raise TypeError("Gene can take only positive integer values")
        if x_cor < 0 or y_cor < 0:
            raise TypeError("Gene can take only positive integer values")

        self.x_cor = x_cor
        self.y_cor = y_cor

    def equals(self, point):
        if not isinstance(point, Gene):
            raise TypeError("Objects of type Gene is allowed only")

        is_same = False
        if self.x_cor == point.x_cor and self.y_cor == point.y_cor:
            is_same = True
        return is_same

    def print(self, message):
        if not isinstance(message, str):
            raise TypeError("Only String message is allowed.")

        print(message, " (", self.x_cor, ", ", self.y_cor, ")")



