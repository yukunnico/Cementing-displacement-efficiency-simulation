const pptxgen = require('pptxgenjs');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'Sisyphus';
pptx.company = 'OhMyOpenCode';
pptx.subject = '固井注替阶段顶替效率实时模拟算法研究';
pptx.title = '固井注替阶段顶替效率实时模拟算法研究';
pptx.lang = 'zh-CN';
pptx.theme = {
  headFontFace: 'Microsoft YaHei',
  bodyFontFace: 'Microsoft YaHei',
  lang: 'zh-CN'
};

const COLORS = {
  navy: '0F2747',
  blue: '1D4E89',
  teal: '2B7A78',
  mint: 'EAF7F6',
  light: 'F7FAFC',
  orange: 'D97706',
  orangeLight: 'FFF7ED',
  green: '1F7A4C',
  greenLight: 'F0FDF4',
  gray: '475569',
  gray2: '64748B',
  border: 'D6E0EA',
  white: 'FFFFFF'
};

function addHeader(slide, title, subtitle) {
  slide.addShape(pptx.ShapeType.rect, {
    x: 0, y: 0, w: 13.333, h: 0.75,
    line: { color: COLORS.navy, transparency: 100 },
    fill: { color: COLORS.navy }
  });
  slide.addText(title, {
    x: 0.55, y: 0.16, w: 8.6, h: 0.26,
    fontFace: 'Microsoft YaHei', fontSize: 24, bold: true,
    color: COLORS.white, margin: 0
  });
  slide.addText(subtitle, {
    x: 9.2, y: 0.18, w: 3.55, h: 0.22,
    fontFace: 'Microsoft YaHei', fontSize: 10,
    color: 'D9E3F0', align: 'right', margin: 0
  });
}

function addPanel(slide, x, y, w, h, title, titleColor = COLORS.blue, fill = COLORS.white) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h,
    rectRadius: 0.08,
    line: { color: COLORS.border, pt: 1.2 },
    fill: { color: fill }
  });
  slide.addText(title, {
    x: x + 0.18, y: y + 0.12, w: w - 0.36, h: 0.22,
    fontFace: 'Microsoft YaHei', fontSize: 16, bold: true,
    color: titleColor, margin: 0
  });
}

function addFormulaBox(slide, x, y, w, h, lines, fillColor) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h,
    rectRadius: 0.06,
    line: { color: COLORS.border, pt: 1 },
    fill: { color: fillColor }
  });
  slide.addText(lines.join('\n'), {
    x: x + 0.16, y: y + 0.12, w: w - 0.32, h: h - 0.2,
    fontFace: 'Consolas', fontSize: 15, bold: true,
    color: COLORS.navy, margin: 0,
    breakLine: false, fit: 'shrink'
  });
}

function addBulletText(slide, x, y, w, items, fontSize = 13) {
  let cursorY = y;
  items.forEach((item) => {
    slide.addText('• ' + item, {
      x, y: cursorY, w, h: 0.24,
      fontFace: 'Microsoft YaHei', fontSize,
      color: COLORS.gray, margin: 0,
      fit: 'shrink'
    });
    cursorY += 0.27;
  });
}

function addFlowBox(slide, x, y, w, h, title, body, fillColor, titleColor = COLORS.navy) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h,
    rectRadius: 0.06,
    line: { color: titleColor, pt: 1.2 },
    fill: { color: fillColor }
  });
  slide.addText(title, {
    x: x + 0.12, y: y + 0.12, w: w - 0.24, h: 0.22,
    fontFace: 'Microsoft YaHei', fontSize: 14, bold: true,
    color: titleColor, margin: 0, align: 'center'
  });
  slide.addText(body, {
    x: x + 0.12, y: y + 0.42, w: w - 0.24, h: h - 0.5,
    fontFace: 'Microsoft YaHei', fontSize: 11,
    color: COLORS.gray, margin: 0, align: 'center', valign: 'mid',
    fit: 'shrink'
  });
}

function addArrowText(slide, x, y, text = '→') {
  slide.addText(text, {
    x, y, w: 0.35, h: 0.22,
    fontFace: 'Microsoft YaHei', fontSize: 20, bold: true,
    color: COLORS.gray2, align: 'center', margin: 0
  });
}

function addVerticalArrow(slide, x, y, h = 0.3) {
  slide.addText('↓', {
    x, y, w: 0.3, h,
    fontFace: 'Microsoft YaHei', fontSize: 20, bold: true,
    color: COLORS.gray2, align: 'center', margin: 0
  });
}

// Slide 1
{
  const slide = pptx.addSlide();
  slide.background = { color: COLORS.light };
  addHeader(slide, '实时顶替效率模拟算法', '依据：固井顶替效率模型说明.md / 流程图.svg / hu102_tail_d2dga_model.py');

  addPanel(slide, 0.45, 0.95, 5.65, 6.0, '左侧：实时效率公式与变量说明', COLORS.blue, COLORS.white);
  addPanel(slide, 6.3, 0.95, 6.55, 6.0, '右侧：模型求解流程图', COLORS.teal, COLORS.white);

  addFormulaBox(slide, 0.7, 1.45, 5.15, 0.75, [
    'cement = x_lead + x_tail'
  ], COLORS.mint);

  addFormulaBox(slide, 0.7, 2.35, 5.15, 1.0, [
    'η(t,y,s) = cement(t,y,s)',
    '           × (1 - wall(t,y,s))'
  ], COLORS.greenLight);

  addFormulaBox(slide, 0.7, 3.55, 5.15, 0.92, [
    'wall(t+dt) = wall(t)',
    '            × exp(-kclean × dt / 150)'
  ], COLORS.orangeLight);

  slide.addText('变量速查', {
    x: 0.72, y: 4.7, w: 1.1, h: 0.2,
    fontFace: 'Microsoft YaHei', fontSize: 14, bold: true,
    color: COLORS.navy, margin: 0
  });

  addBulletText(slide, 0.78, 4.98, 4.9, [
    'x_lead / x_tail：领浆与尾浆体积分数',
    'cement：某一时刻、某一位置的总水泥浆体积分数',
    'wall：壁面泥饼残余率；1表示未清洗，0表示完全清洗',
    'kclean：综合清洗速率，dt：时间步长',
    'η(t,y,s)：该位置、该时刻的实时顶替效率'
  ], 12.5);

  addFlowBox(slide, 6.55, 1.55, 1.35, 0.9, '1 输入参数', '环空几何\n排量\n注入时序', 'F8FAFC', COLORS.blue);
  addFlowBox(slide, 8.15, 1.55, 1.55, 0.9, '2 数值模拟推进', '对流 + 扩散\n逐时间步更新', 'EFF6FF', COLORS.blue);
  addFlowBox(slide, 10.0, 1.55, 1.6, 0.9, '3 水泥浆浓度', 'cement =\nx_lead + x_tail', 'EFF6FF', COLORS.teal);
  addFlowBox(slide, 10.25, 3.05, 1.45, 0.95, '4 壁面清洗', 'wall 衰减\nexp(-kclean·dt/150)', 'FFF7ED', COLORS.orange);
  addFlowBox(slide, 7.75, 3.0, 2.0, 1.1, '5 实时顶替效率', 'η(t,y,s) = cement × (1 - wall)', 'F0FDF4', COLORS.green);

  addArrowText(slide, 7.92, 1.82);
  addArrowText(slide, 9.77, 1.82);
  addVerticalArrow(slide, 10.82, 2.48);
  addArrowText(slide, 9.86, 3.42, '←');

  slide.addText('实时算法要点：先更新局部水泥浆浓度与壁面泥饼残余率，再在每个时空网格点计算 η(t,y,s)。', {
    x: 6.6, y: 5.15, w: 6.0, h: 0.35,
    fontFace: 'Microsoft YaHei', fontSize: 12.5,
    color: COLORS.gray, margin: 0
  });

  slide.addText('代码依据：hu102_tail_d2dga_model.py 第705–720行', {
    x: 6.6, y: 5.55, w: 4.0, h: 0.18,
    fontFace: 'Microsoft YaHei', fontSize: 10.5,
    color: COLORS.gray2, margin: 0
  });
}

// Slide 2
{
  const slide = pptx.addSlide();
  slide.background = { color: COLORS.light };
  addHeader(slide, '最终顶替效率模拟', '依据：固井顶替效率模型说明.md / 流程图.svg / hu102_tail_d2dga_model.py');

  addPanel(slide, 0.45, 0.95, 5.65, 6.0, '左侧：最终效率公式与物理意义', COLORS.blue, COLORS.white);
  addPanel(slide, 6.3, 0.95, 6.55, 6.0, '右侧：空间加权积分流程', COLORS.teal, COLORS.white);

  addFormulaBox(slide, 0.7, 1.45, 5.15, 1.08, [
    'η_final = ∬(b × η) dy ds',
    '          ───────────────',
    '              ∬b dy ds'
  ], COLORS.greenLight);

  addFormulaBox(slide, 0.7, 2.72, 5.15, 1.02, [
    'η_interval = ∬(b × η × mask) dy ds',
    '             ─────────────────────',
    '               ∬(b × mask) dy ds'
  ], COLORS.mint);

  slide.addText('物理意义', {
    x: 0.72, y: 3.98, w: 1.1, h: 0.2,
    fontFace: 'Microsoft YaHei', fontSize: 14, bold: true,
    color: COLORS.navy, margin: 0
  });
  addBulletText(slide, 0.78, 4.26, 4.92, [
    '将全空间内每个位置的实时顶替效率 η 按局部环空间隙 b 加权平均',
    '全井段、目标井段和 CBL 井段的区别仅在于积分时采用的 mask 不同',
    '若局部间隙更大，则该位置对最终效率的贡献更高'
  ], 12.5);

  slide.addText('变量速查', {
    x: 0.72, y: 5.25, w: 1.1, h: 0.2,
    fontFace: 'Microsoft YaHei', fontSize: 14, bold: true,
    color: COLORS.navy, margin: 0
  });
  addBulletText(slide, 0.78, 5.53, 4.9, [
    'b：局部环空间隙',
    'η：局部实时顶替效率',
    'η_final：最终顶替效率',
    'mask：目标井段或 CBL 井段掩码'
  ], 12.5);

  addFlowBox(slide, 6.7, 1.6, 2.15, 1.1, '1 实时顶替效率场', 'η(t,y,s) 已由前一阶段求得', 'F0FDF4', COLORS.green);
  addFlowBox(slide, 9.25, 1.5, 2.3, 1.28, '2 加权积分', '对整个环空或目标井段\n计算 ∬(b·η) dy ds\n并用 ∬b dy ds 归一化', 'EFF6FF', COLORS.blue);
  addFlowBox(slide, 9.2, 3.4, 2.4, 1.15, '3 最终顶替效率', '输出 η_final 或 η_interval', 'FFF7ED', COLORS.orange);

  addArrowText(slide, 8.92, 1.98);
  addVerticalArrow(slide, 10.2, 2.95);

  slide.addText('不同井段仅 mask 不同，公式结构完全一致：\n全井段 η_final、目标井段 η_target、CBL井段 η_CBL 可统一到同一积分框架。', {
    x: 6.7, y: 5.15, w: 5.6, h: 0.5,
    fontFace: 'Microsoft YaHei', fontSize: 12.5,
    color: COLORS.gray, margin: 0
  });

  slide.addText('代码依据：hu102_tail_d2dga_model.py 第722–728行；说明文档第40–62行', {
    x: 6.7, y: 5.75, w: 5.4, h: 0.18,
    fontFace: 'Microsoft YaHei', fontSize: 10.5,
    color: COLORS.gray2, margin: 0
  });
}

async function main() {
  const out = 'D:/users/desktop/research/海能发项目相关文件/相关PPT/固井注替阶段顶替效率模拟算法研究.pptx';
  await pptx.writeFile({ fileName: out });
  console.log(`PPT generated: ${out}`);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
