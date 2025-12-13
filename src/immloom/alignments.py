import re
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

# ПАРСЕР: принимает либо путь к файлу, либо текст выравнивания
def parse_clustal(text: str):
    """
    Возвращает dict: {seq_id: full_sequence}
    Игнорирует заголовки и номера в конце строки.
    """
    seqs = defaultdict(list)
    # регулярное выражение: id <whitespace> sequence_fragment [opt number]
    pattern = re.compile(r'^(\S+)\s+([A-Za-z\-]+)\s*(?:\d+)?$')
    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        if not line:
            continue
        if line.startswith('CLUSTAL') or line.startswith('#') or line.startswith(' '):
            # пропускаем заголовок и общие строки
            continue
        m = pattern.match(line)
        if m:
            sid = m.group(1)
            frag = m.group(2)
            seqs[sid].append(frag)
    # Собираем полные последовательности
    full = {sid: ''.join(frags) for sid, frags in seqs.items()}
    return full

# ВЫЧИСЛЕНИЕ позиций first/last (1-based). Если нет ни одной ненулевой — ставим first=0,last=0
def compute_intervals(seqs_dict):
    rows = []
    all_lengths = {len(s) for s in seqs_dict.values()}
    if len(all_lengths) > 1:
        # предупреждение (можно выровнять до максимума, но лучше знать).
        print(f"Warning: different alignment fragment lengths found (distinct lengths: {sorted(all_lengths)})")
    align_len = max(all_lengths) if all_lengths else 0
    for sid, seq in seqs_dict.items():
        first = 0
        last = 0
        if seq:
            for i, ch in enumerate(seq):
                if ch != '-':
                    first = i + 1
                    break
            if first != 0:
                for i, ch in enumerate(reversed(seq)):
                    if ch != '-':
                        last = align_len - i
                        break
        num_gaps = seq.count('-')
        rows.append({
            'seq_id': sid,
            'first': first,
            'last': last,
            'segment_length': (last - first + 1) if (first and last) else 0,
            'alignment_length': align_len,
            'num_gaps': num_gaps,
        })
    df = pd.DataFrame(rows).sort_values(['alignment_length','seq_id'], ascending=[False, True]).reset_index(drop=True)
    return df

# ВИЗУАЛИЗАЦИЯ: один отрезок (first..last) на строку
def plot_intervals(df, title='Segments (first..last) on alignment'):
    # сортируем по порядку для приятного отображения
    df = df.copy()
    df['plot_order'] = range(len(df), 0, -1)  # so top = first row
    fig_h = max(4, 0.25 * len(df))
    fig, ax = plt.subplots(figsize=(12, fig_h))
    yticks = []
    yticklabels = []
    for i, row in df.iterrows():
        y = i
        yticks.append(y)
        yticklabels.append(row['seq_id'])
        if row['first'] and row['last'] and row['last'] >= row['first']:
            # broken_barh expects (xmin, width)
            ax.broken_barh([(row['first'], row['last'] - row['first'] + 1)], (y - 0.4, 0.8), facecolors=('tab:blue',))
        else:
            # пустые последовательности — рисуем пунктирную линию/mark
            ax.plot([], [])  # noop to keep consistent
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=8)
    ax.set_xlabel('Alignment position (1-based)')
    ax.set_title(title)
    ax.set_xlim(1, int(df['alignment_length'].max() if not df.empty else 1))
    ax.invert_yaxis()
    ax.grid(axis='x', linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.show()

# Удобные обёртки:
def alignments_plot(path):
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read()
    seqs = parse_clustal(text)
    df = compute_intervals(seqs)
    # df.to_csv('segments.csv', index=False)
    # print("Saved segments.csv (columns: seq_id, first, last, segment_length, alignment_length, num_gaps).")
    display(df)  # если в jupyter — покажет табличку
    plot_intervals(df)
    return df

def from_text_and_plot(alignment_text):
    seqs = parse_clustal(alignment_text)
    df = compute_intervals(seqs)
    # df.to_csv('segments.csv', index=False)
    # print("Saved segments.csv (columns: seq_id, first, last, segment_length, alignment_length, num_gaps).")
    display(df)
    plot_intervals(df)
    return df