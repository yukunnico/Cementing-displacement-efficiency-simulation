const pptxgen = require('pptxgenjs');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'Sisyphus';
pptx.company = 'OhMyOpenCode';
pptx.subject = '固井注替阶段顶替效率实时模拟算法研究（论文补充版）';
pptx.title = '固井注替阶段顶替效率实时模拟算法研究（论文补充版）';
pptx.lang = 'zh-CN';
pptx.theme = {
  headFontFace: 'Microsoft YaHei',
  bodyFontFace: 'Microsoft YaHei',
  lang: 'zh-CN'
};

const C = {
  navy: '0F2747', blue: '1D4E89', teal: '2B7A78', green: '1F7A4C', orange: 'D97706',
  mint: 'EAF7F6', light: 'F7FAFC', greenLight: 'F0FDF4', orangeLight: 'FFF7ED', pinkLight: 'FDF2F8',
  white: 'FFFFFF', gray: '475569', gray2: '64748B', border: 'D6E0EA'
};

function header(slide, title, sub) {
  slide.addShape(pptx.ShapeType.rect, { x:0, y:0, w:13.333, h:0.75, line:{color:C.navy, transparency:100}, fill:{color:C.navy} });
  slide.addText(title, { x:0.55, y:0.16, w:8.8, h:0.3, fontFace:'Microsoft YaHei', fontSize:24, bold:true, color:C.white, margin:0 });
  slide.addText(sub, { x:8.9, y:0.18, w:3.9, h:0.22, fontFace:'Microsoft YaHei', fontSize:10, color:'D9E3F0', align:'right', margin:0, fit:'shrink' });
}

function panel(slide, x, y, w, h, title, color, fill=C.white) {
  slide.addShape(pptx.ShapeType.roundRect, { x, y, w, h, rectRadius:0.08, line:{color:C.border, pt:1.2}, fill:{color:fill} });
  slide.addText(title, { x:x+0.18, y:y+0.12, w:w-0.36, h:0.22, fontFace:'Microsoft YaHei', fontSize:16, bold:true, color, margin:0 });
}

function fbox(slide, x, y, w, h, lines, fill) {
  slide.addShape(pptx.ShapeType.roundRect, { x, y, w, h, rectRadius:0.05, line:{color:C.border, pt:1}, fill:{color:fill} });
  slide.addText(lines.join('\n'), { x:x+0.15, y:y+0.1, w:w-0.3, h:h-0.18, fontFace:'Consolas', fontSize:13.5, bold:true, color:C.navy, margin:0, fit:'shrink' });
}

function bullets(slide, x, y, w, items, fs=11.8, step=0.24) {
  let cy=y;
  items.forEach(item => {
    slide.addText('• ' + item, { x, y:cy, w, h:0.22, fontFace:'Microsoft YaHei', fontSize:fs, color:C.gray, margin:0, fit:'shrink' });
    cy += step;
  });
}

function flow(slide, x, y, w, h, title, body, fill, color=C.navy) {
  slide.addShape(pptx.ShapeType.roundRect, { x, y, w, h, rectRadius:0.05, line:{color, pt:1.2}, fill:{color:fill} });
  slide.addText(title, { x:x+0.12, y:y+0.11, w:w-0.24, h:0.2, fontFace:'Microsoft YaHei', fontSize:13.5, bold:true, color, margin:0, align:'center' });
  slide.addText(body, { x:x+0.12, y:y+0.37, w:w-0.24, h:h-0.46, fontFace:'Microsoft YaHei', fontSize:10.8, color:C.gray, margin:0, align:'center', valign:'mid', fit:'shrink' });
}

function arrow(slide, x, y, text='→') {
  slide.addText(text, { x, y, w:0.35, h:0.22, fontFace:'Microsoft YaHei', fontSize:20, bold:true, color:C.gray2, align:'center', margin:0 });
}

function down(slide, x, y) {
  slide.addText('↓', { x, y, w:0.3, h:0.25, fontFace:'Microsoft YaHei', fontSize:20, bold:true, color:C.gray2, align:'center', margin:0 });
}

// Slide 1
{
  const s = pptx.addSlide();
  s.background = { color: C.light };
  header(s, '实时顶替效率模拟算法（论文基础 + 工程扩展）', '论文：Zhang & Frigaard (2022)；代码：hu102_tail_d2dga_model.py');
  panel(s, 0.45, 0.95, 5.8, 6.0, '左侧：论文二维基础方程 + 代码实时效率公式', C.blue);
  panel(s, 6.35, 0.95, 6.45, 6.0, '右侧：模型求解流程图', C.teal);

  fbox(s, 0.7, 1.35, 5.3, 1.02, [
    '∂/∂t [H r_a c̄] + ∇_a·q = 0',
    '∇_a·[S + b] = 0'
  ], C.mint);
  s.addText('论文 D2DGA 二维基础：体积分数输运 + 流函数椭圆方程', { x:0.74, y:2.43, w:5.1, h:0.18, fontFace:'Microsoft YaHei', fontSize:11.2, color:C.gray2, margin:0 });

  fbox(s, 0.7, 2.75, 5.3, 0.78, [
    'cement = x_lead + x_tail'
  ], C.greenLight);
  fbox(s, 0.7, 3.68, 5.3, 0.94, [
    'wall(t+dt) = wall(t) × exp(-kclean × dt / 150)'
  ], C.orangeLight);
  fbox(s, 0.7, 4.78, 5.3, 0.95, [
    'η(t,y,s) = cement(t,y,s) × (1 - wall(t,y,s))'
  ], C.pinkLight);

  bullets(s, 0.78, 5.82, 5.05, [
    '论文提供二维顶替与分散方程；代码在此基础上增加壁面泥饼清除与实时效率定义。',
    '代码依据：第705–720行，先求 cement、再更新 wall、最后计算 eff。'
  ], 11.4, 0.23);

  flow(s, 6.62, 1.45, 1.25, 0.86, '1 输入参数', '环空几何\n排量\n注入时序', 'F8FAFC', C.blue);
  flow(s, 8.08, 1.45, 1.45, 0.86, '2 数值模拟推进', '对流 + 扩散\n逐时间步更新', 'EFF6FF', C.blue);
  flow(s, 9.78, 1.45, 1.55, 0.86, '3 水泥浆浓度', 'cement =\nx_lead + x_tail', 'EFF6FF', C.teal);
  flow(s, 10.03, 2.95, 1.45, 0.92, '4 壁面清洗', 'wall 衰减\nexp(-kclean·dt/150)', 'FFF7ED', C.orange);
  flow(s, 7.75, 2.95, 2.0, 1.05, '5 实时顶替效率', 'η(t,y,s) = cement × (1 - wall)', 'F0FDF4', C.green);
  arrow(s, 7.88, 1.72); arrow(s, 9.58, 1.72); down(s, 10.58, 2.38); arrow(s, 9.68, 3.31, '←');
  s.addText('右图表达的是当前工程化求解流程；其底层二维输运思想与论文 D2DGA 一致，但局部效率定义来自代码扩展。', {
    x:6.58, y:5.18, w:6.0, h:0.45, fontFace:'Microsoft YaHei', fontSize:12, color:C.gray, margin:0
  });
}

// Slide 2
{
  const s = pptx.addSlide();
  s.background = { color: C.light };
  header(s, '最终顶替效率模拟（论文定义 vs 代码定义）', '重点：论文 η_E 与代码 η_final 的差异和联系');
  panel(s, 0.45, 0.95, 5.8, 6.0, '左侧：论文与代码的最终效率公式', C.blue);
  panel(s, 6.35, 0.95, 6.45, 6.0, '右侧：最终效率形成过程', C.teal);

  fbox(s, 0.7, 1.32, 5.3, 0.9, [
    '论文：η_E = V_displaced(t*) / V_annulus',
    '      t* = 1.2 × L_annulus / w_0'
  ], C.mint);
  s.addText('论文中的位移效率更偏“总体置换体积分数”定义。', { x:0.74, y:2.29, w:5.0, h:0.18, fontFace:'Microsoft YaHei', fontSize:11.2, color:C.gray2, margin:0 });

  fbox(s, 0.7, 2.68, 5.3, 1.0, [
    '代码：η_final = ∬(b × η) dy ds / ∬b dy ds'
  ], C.greenLight);
  fbox(s, 0.7, 3.86, 5.3, 0.96, [
    'η_interval = ∬(b × η × mask) dy ds',
    '             / ∬(b × mask) dy ds'
  ], C.orangeLight);

  bullets(s, 0.78, 5.0, 5.1, [
    '论文 η_E：看“置换了多少体积”。',
    '代码 η_final：先算局部有效顶替效率 η，再按环空间隙 b 做空间加权平均。',
    '因此代码的最终效率比论文原始 η_E 更适合表达局部清洗与井段评价。'
  ], 11.6, 0.25);

  flow(s, 6.72, 1.55, 2.0, 1.0, '1 论文二维浓度场', '由 2DGA / D2DGA\n得到 c̄(φ,ξ,t)', 'EFF6FF', C.blue);
  flow(s, 9.0, 1.55, 2.15, 1.0, '2 代码局部效率场', 'η = cement × (1 - wall)', 'F0FDF4', C.green);
  flow(s, 8.9, 3.2, 2.35, 1.12, '3 空间加权积分', '对全井段或目标井段\n计算 ∬(b·η) / ∬b', 'FFF7ED', C.orange);
  flow(s, 8.82, 4.75, 2.5, 1.0, '4 最终效率输出', 'η_final / η_target / η_CBL', 'FDF2F8', C.teal);
  arrow(s, 8.76, 1.93); down(s, 9.98, 2.68); down(s, 9.98, 4.3);
  s.addText('右图说明：代码并未直接采用论文 η_E 作为最终指标，而是在论文二维模型思想上叠加了局部有效效率与井段加权积分。', {
    x:6.62, y:6.0, w:5.95, h:0.38, fontFace:'Microsoft YaHei', fontSize:12, color:C.gray, margin:0
  });
}

async function main() {
  const out = 'D:/users/desktop/research/海能发项目相关文件/相关PPT/固井注替阶段顶替效率模拟算法研究-论文补充版.pptx';
  await pptx.writeFile({ fileName: out });
  console.log(`PPT generated: ${out}`);
}

main().catch(err => { console.error(err); process.exit(1); });
