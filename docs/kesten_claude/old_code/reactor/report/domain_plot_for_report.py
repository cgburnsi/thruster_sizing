import matplotlib.pyplot as plt

# Define regions
x_total = []
y_total = []

ADIV = 100
BDIV = 100
A = 1.0
B = 0.3

# Region 1: 0 to B
x1 = [i * B / ADIV for i in range(ADIV + 1)]
y1 = [0.5 for _ in x1]

# Region 2: B to A
x2 = [B + i * (A - B) / BDIV for i in range(1, BDIV + 1)]
y2 = [1.0 for _ in x2]

# Combine
x_total = x1 + x2
y_total = y1 + y2

# Plot
fig, ax = plt.subplots(figsize=(10, 2))
ax.plot(x1, y1, 'b.-', label="Region 1: 0 to B (DR)")
ax.plot(x2, y2, 'r.-', label="Region 2: B to A (BDR)")
ax.axvline(B, color='k', linestyle='--', label='Transition at B')
ax.set_yticks([])
ax.set_xlabel("Domain Position (x)")
ax.set_title("Domain Discretization in SLOPE Subroutine")
ax.legend(loc='upper center', ncol=3)
ax.grid(True)

plt.tight_layout()
plt.show()
